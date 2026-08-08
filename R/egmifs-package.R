#' egmifs: Extended generalized monotone incremental forward stagewise regression
#'
#' A modular R/Rcpp/C++ framework for extended generalized monotone incremental
#' forward stagewise regression. The core path engine is separated from family,
#' link, fused family-link, and information-criterion modules. Built-in modules
#' provide NB2 and Poisson models, while custom differentiable modules can reuse
#' the same stagewise path construction.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom nloptr nloptr
#' @useDynLib egmifs, .registration = TRUE
## usethis namespace: end
NULL
