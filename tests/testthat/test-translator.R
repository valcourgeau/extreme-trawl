test_that("standard to transform - positive xi", {
  params <- c(.5, 1, .5)
  parametrisation <- "standard"
  target <- "transform"
  translation <- parametrisation_translator(
    params = params,
    parametrisation = parametrisation,
    target = target
  )
  testthat::expect_equal(
    translation,
    c(1., 1. + params[3], params[3])
  )
})

test_that("standard to transform - negative xi", {
  params <- c(-.5, 1, .5)
  parametrisation <- "standard"
  target <- "transform"
  translation <- parametrisation_translator(
    params = params,
    parametrisation = parametrisation,
    target = target
  )
  testthat::expect_equal(
    translation,
    c(1., 1. + params[3], params[3])
  )
})


test_that("standard to noven - positive xi", {
  params <- c(.5, 1, .5)
  parametrisation <- "standard"
  target <- "noven"
  translation <- parametrisation_translator(
    params = params,
    parametrisation = parametrisation,
    target = target
  )
  testthat::expect_equal(
    translation,
    c(1. / params[1], params[2] / params[1] - params[3], params[3])
  )
})

test_that("standard to noven - negative xi", {
  params <- c(-.5, 1, .5)
  parametrisation <- "standard"
  target <- "noven"
  translation <- parametrisation_translator(
    params = params,
    parametrisation = parametrisation,
    target = target
  )
  testthat::expect_equal(
    translation,
    c(1. / params[1], params[2] / abs(params[1]) - params[3], params[3])
  )
})

test_that("noven to standard - positive xi", {
  params <- c(15., 5., .5)
  parametrisation <- "noven"
  target <- "standard"
  translation <- parametrisation_translator(
    params = params,
    parametrisation = parametrisation,
    target = target
  )
  testthat::expect_equal(
    translation,
    c(1. / params[1], (params[2] + params[3]) / params[1], params[3])
  )
})

test_that("noven to standard - negative xi", {
  params <- c(-15., 5., .5)
  parametrisation <- "noven"
  target <- "standard"
  translation <- parametrisation_translator(
    params = params,
    parametrisation = parametrisation,
    target = target
  )
  testthat::expect_equal(
    translation,
    c(1. / params[1], (params[2] + params[3]) / abs(params[1]), params[3])
  )
})

test_that("same to same - negative xi", {
  params <- c(-15., 5., .5)
  parametrisation <- "standard"
  target <- "standard"
  translation <- parametrisation_translator(
    params = params, parametrisation = parametrisation, target = target
  )
  testthat::expect_equal(translation, translation)
})

test_that("noven to transform - negative xi", {
  params <- c(-15., 5., .5)
  parametrisation <- "noven"
  target <- "transform"
  translation <- parametrisation_translator(
    params = params, parametrisation = parametrisation, target = target
  )
  testthat::expect_equal(
    translation,
    c(1., (1.0 + params[3]) / abs(1.0), params[3])
  )
})

test_that("wrong parametrisation", {
  params <- c(-15., 5., .5)
  parametrisation <- "transform"
  target <- "standard"
  testthat::expect_error(parametrisation_translator(
    params = params, parametrisation = parametrisation, target = target
  ), regexp = "is not TRUE")
})


test_that("wrong parametrisation error", {
  params <- c(-15., 5., .5)
  parametrisation <- "transform"
  target <- "tmp"
  testthat::expect_error(parametrisation_translator(
    params = params, parametrisation = parametrisation, target = target
  ), regexp = "is not TRUE")
})
