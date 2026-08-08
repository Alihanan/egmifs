.egmifs.call.argument <- function(object, name) {
  model.call <- object$call

  if (is.null(model.call) || is.null(model.call[[name]])) {
    return(NULL)
  }

  model.call[[name]]
}


.egmifs.expression.function.name <- function(expression) {
  if (is.null(expression)) {
    return(NA_character_)
  }

  if (is.character(expression) && length(expression) == 1L) {
    return(as.character(expression))
  }

  if (is.name(expression)) {
    return(as.character(expression))
  }

  if (!is.call(expression) || length(expression) == 0L) {
    return(NA_character_)
  }

  head <- expression[[1L]]

  if (is.name(head)) {
    return(as.character(head))
  }

  if (
      is.call(head) &&
      length(head) == 3L &&
      is.name(head[[1L]]) &&
      as.character(head[[1L]]) %in% c("::", ":::")
  ) {
    return(
      paste0(
        as.character(head[[2L]]),
        as.character(head[[1L]]),
        as.character(head[[3L]])
      )
    )
  }

  NA_character_
}


.egmifs.resolve.function <- function(name, envir) {
  if (
      length(name) != 1L ||
      is.na(name) ||
      !nzchar(name)
  ) {
    return(NULL)
  }

  namespace.match <- regexec(
    "^([^:]+)(:::{0,1})([^:]+)$",
    name
  )
  namespace.parts <- regmatches(name, namespace.match)[[1L]]

  if (length(namespace.parts) == 4L) {
    package <- namespace.parts[[2L]]
    operator <- namespace.parts[[3L]]
    function.name <- namespace.parts[[4L]]

    return(
      tryCatch(
        if (identical(operator, "::")) {
          getExportedValue(package, function.name)
        } else {
          get(function.name, envir = asNamespace(package), inherits = FALSE)
        },
        error = function(...) NULL
      )
    )
  }

  candidates <- list(
    envir,
    asNamespace("egmifs"),
    .GlobalEnv
  )

  for (candidate in candidates) {
    value <- get0(
      name,
      envir = candidate,
      mode = "function",
      inherits = TRUE,
      ifnotfound = NULL
    )

    if (is.function(value)) {
      return(value)
    }
  }

  NULL
}


.egmifs.call.metadata <- function(expression, envir) {
  function.name <- .egmifs.expression.function.name(expression)

  resolved.function <- if (is.character(expression)) {
    NULL
  } else {
    .egmifs.resolve.function(function.name, envir)
  }

  list(
    call = expression,
    function_name = function.name,
    constructor = resolved.function
  )
}


.egmifs.named.plugin.metadata <- function(value, fallback.name) {
  if (is.list(value) && length(value) > 0L) {
    value.names <- names(value)
    name <- if (
        !is.null(value.names) &&
        length(value.names) > 0L &&
        !is.na(value.names[[1L]]) &&
        nzchar(value.names[[1L]])
    ) {
      value.names[[1L]]
    } else {
      fallback.name
    }

    details <- value[[1L]]
    parameter.count <- if (
        is.list(details) &&
        !is.null(details$parameter_count)
    ) {
      as.integer(details$parameter_count[[1L]])
    } else {
      NA_integer_
    }

    return(
      list(
        name = name,
        parameter_count = parameter.count,
        metadata = details
      )
    )
  }

  list(
    name = fallback.name,
    parameter_count = NA_integer_,
    metadata = value
  )
}


.egmifs.component.metadata <- function(
    object,
    component,
    envir
) {
  component <- match.arg(component, c("family", "link"))
  input <- .egmifs.input.view(object)

  family.value <- if (is.list(input)) input$family else input
  link.value <- if (is.list(input)) input$link_func else input
  family.link.supplied <-
    is.list(input) && isTRUE(input$family_link_supplied)

  if (identical(component, "family")) {
    plugin <- .egmifs.named.plugin.metadata(
      family.value,
      .egmifs.family.label(object)
    )
  } else {
    plugin <- .egmifs.named.plugin.metadata(
      link.value,
      .egmifs.link.label(object)
    )
  }

  family.link.call <- .egmifs.call.argument(object, "family.link")
  component.call <- .egmifs.call.argument(object, component)

  specification <- if (!is.null(family.link.call)) {
    family.link.call
  } else {
    component.call
  }

  source <- if (!is.null(family.link.call)) {
    "family.link"
  } else if (!is.null(component.call)) {
    component
  } else if (family.link.supplied) {
    "automatic family.link"
  } else {
    "default"
  }

  call.metadata <- .egmifs.call.metadata(specification, envir)

  out <- c(
    list(
      name = plugin$name,
      parameter_count = plugin$parameter_count,
      mode = .egmifs.interface.label(object),
      source = source,
      family_link_supplied = family.link.supplied,
      metadata = plugin$metadata
    ),
    call.metadata
  )

  class(out) <- c(
    paste0("egmifs_", component),
    "egmifs_component",
    "list"
  )
  out
}


.egmifs.metadata.collection <- function(
    values,
    labels,
    drop
) {
  names(values) <- labels

  if (length(values) == 1L && isTRUE(drop)) {
    return(values[[1L]])
  }

  values
}


#' Extract the fitted family specification
#'
#' The standard [stats::family()] generic returns compact egmifs family
#' metadata for a fitted model. The result includes the family name, parameter
#' count, separate-versus-fused interface mode, original call specification,
#' constructor name, and the constructor function itself when it can
#' still be resolved in the calling environment or a namespace.
#'
#' @param object An egmifs fit, multi-alpha fit, or tuning grid.
#' @param ... Additional arguments passed to the class method.
#'
#' @return A named metadata list for one fit, or a named list of metadata lists
#'   when multiple alpha or alpha-by-eta fits are selected.
#' @export
family.egmifs <- function(object, ...) {
  .egmifs.component.metadata(
    object,
    "family",
    envir = parent.frame()
  )
}


#' @rdname family.egmifs
#' @param alpha Optional stored alpha value or values, or their labels.
#' @param drop Logical. Drop a one-fit result to its metadata list.
#' @export
family.egmifs_multi_alpha <- function(
    object,
    alpha = NULL,
    drop = TRUE,
    ...
) {
  indices <- .egmifs.multi.alpha.select(object, alpha = alpha)
  values <- lapply(
    indices,
    function(index) stats::family(object$fits[[index]], ...)
  )

  .egmifs.metadata.collection(
    values,
    paste0("alpha=", object$alpha_labels[indices]),
    drop
  )
}


#' @rdname family.egmifs
#' @param eta Optional stored prior-strength value or values, or their labels.
#' @export
family.egmifs_tuning_grid <- function(
    object,
    alpha = NULL,
    eta = NULL,
    drop = TRUE,
    ...
) {
  indices <- .egmifs.tuning.grid.select(
    object,
    alpha = alpha,
    eta = eta
  )
  values <- lapply(
    indices,
    function(index) stats::family(object$fits[[index]], ...)
  )

  .egmifs.metadata.collection(
    values,
    names(object$fits)[indices],
    drop
  )
}


#' Extract the fitted link specification
#'
#' `link()` reports the link name, parameter count, interface mode, original
#' call specification, constructor name, and the constructor function
#' when it can still be resolved.
#'
#' @inheritParams family.egmifs
#'
#' @return A named metadata list for one fit, or a named list of metadata lists
#'   when multiple fits are selected.
#' @export
link <- function(object, ...) {
  UseMethod("link")
}


#' @rdname link
#' @export
link.egmifs <- function(object, ...) {
  .egmifs.component.metadata(
    object,
    "link",
    envir = parent.frame()
  )
}


#' @rdname link
#' @inheritParams family.egmifs_multi_alpha
#' @export
link.egmifs_multi_alpha <- function(
    object,
    alpha = NULL,
    drop = TRUE,
    ...
) {
  indices <- .egmifs.multi.alpha.select(object, alpha = alpha)
  values <- lapply(
    indices,
    function(index) link(object$fits[[index]], ...)
  )

  .egmifs.metadata.collection(
    values,
    paste0("alpha=", object$alpha_labels[indices]),
    drop
  )
}


#' @rdname link
#' @inheritParams family.egmifs_tuning_grid
#' @export
link.egmifs_tuning_grid <- function(
    object,
    alpha = NULL,
    eta = NULL,
    drop = TRUE,
    ...
) {
  indices <- .egmifs.tuning.grid.select(
    object,
    alpha = alpha,
    eta = eta
  )
  values <- lapply(
    indices,
    function(index) link(object$fits[[index]], ...)
  )

  .egmifs.metadata.collection(
    values,
    names(object$fits)[indices],
    drop
  )
}


.egmifs.criteria.call.specifications <- function(object) {
  expression <- .egmifs.call.argument(object, "criteria")

  if (is.null(expression)) {
    return(list(source_call = NULL, specifications = list()))
  }

  if (
      is.call(expression) &&
      length(expression) >= 1L &&
      is.name(expression[[1L]]) &&
      identical(as.character(expression[[1L]]), "list")
  ) {
    values <- as.list(expression)[-1L]
    return(
      list(
        source_call = expression,
        specifications = values
      )
    )
  }

  list(
    source_call = expression,
    specifications = list()
  )
}


.egmifs.criteria.metadata <- function(object, envir) {
  input <- .egmifs.input.view(object)
  criterion.names <- if (is.list(input)) {
    as.character(input$criteria)
  } else {
    character()
  }

  if (length(criterion.names) == 0L) {
    criterion.names <- names(.egmifs.best.criteria(object))
  }

  call.specifications <- .egmifs.criteria.call.specifications(object)
  specifications <- call.specifications$specifications
  specification.names <- names(specifications)

  entries <- lapply(
    seq_along(criterion.names),
    function(index) {
      specification <- if (index <= length(specifications)) {
        specifications[[index]]
      } else {
        NULL
      }

      call.metadata <- .egmifs.call.metadata(specification, envir)
      label <- if (
          !is.null(specification.names) &&
          index <= length(specification.names) &&
          !is.na(specification.names[[index]]) &&
          nzchar(specification.names[[index]])
      ) {
        specification.names[[index]]
      } else {
        criterion.names[[index]]
      }

      c(
        list(
          name = criterion.names[[index]],
          label = label
        ),
        call.metadata
      )
    }
  )
  names(entries) <- criterion.names

  out <- list(
    criteria = entries,
    source_call = call.specifications$source_call,
    selected = .egmifs.criterion.table(object)
  )
  class(out) <- c("egmifs_criteria", "list")
  out
}


#' Extract fitted criterion specifications
#'
#' `criteria()` returns the criteria used for fitting, their call expressions,
#' constructor names, resolvable constructor functions, and the table
#' of criterion-selected states.
#'
#' @inheritParams family.egmifs
#'
#' @return A metadata list for one fit, or a named list of metadata lists when
#'   multiple fits are selected.
#' @export
criteria <- function(object, ...) {
  UseMethod("criteria")
}


#' @rdname criteria
#' @export
criteria.egmifs <- function(object, ...) {
  .egmifs.criteria.metadata(object, envir = parent.frame())
}


#' @rdname criteria
#' @inheritParams family.egmifs_multi_alpha
#' @export
criteria.egmifs_multi_alpha <- function(
    object,
    alpha = NULL,
    drop = TRUE,
    ...
) {
  indices <- .egmifs.multi.alpha.select(object, alpha = alpha)
  values <- lapply(
    indices,
    function(index) criteria(object$fits[[index]], ...)
  )

  .egmifs.metadata.collection(
    values,
    paste0("alpha=", object$alpha_labels[indices]),
    drop
  )
}


#' @rdname criteria
#' @inheritParams family.egmifs_tuning_grid
#' @export
criteria.egmifs_tuning_grid <- function(
    object,
    alpha = NULL,
    eta = NULL,
    drop = TRUE,
    ...
) {
  indices <- .egmifs.tuning.grid.select(
    object,
    alpha = alpha,
    eta = eta
  )
  values <- lapply(
    indices,
    function(index) criteria(object$fits[[index]], ...)
  )

  .egmifs.metadata.collection(
    values,
    names(object$fits)[indices],
    drop
  )
}
