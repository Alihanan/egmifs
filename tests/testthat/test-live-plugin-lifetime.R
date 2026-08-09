test_that("live plugin pointers survive a second constructor cache miss", {
  skip_on_cran()

  clear.plugin.cache()
  on.exit(
    {
      clear.plugin.cache()
      invisible(gc(verbose = FALSE))
    },
    add = TRUE
  )

  first <- plugin.family.link.nb2.log.live(
    dispersion.initial = 0.5,
    cache = TRUE
  )

  # The C++ source and factory are identical, but the constructor arguments are
  # different, so this is a miss in egmifs' object cache. The first pointer must
  # remain valid when sourceCpp resolves the already-compiled source for the
  # second object.
  second <- plugin.family.link.nb2.log.live(
    dispersion.initial = 1.0,
    cache = TRUE
  )

  first_info <- egmifs:::inspect_family_link_plugin(first)
  second_info <- egmifs:::inspect_family_link_plugin(second)

  expect_identical(first_info$family_name, "NB2.live")
  expect_identical(first_info$link_name, "Log.live")
  expect_identical(second_info$family_name, "NB2.live")
  expect_identical(second_info$link_name, "Log.live")
})
