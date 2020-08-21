test_that("standard to transform - positive xi", {
  params <- c(.5, 1, .5)
  parametrisation <- "standard"
  target <- "transform"
  translation <- ParametrisationTranslator(
    params = params,
    parametrisation = parametrisation,
    target = target)
  expect_equal(translation, c(1., 1. + params[3], params[3]))
})

test_that("standard to transform - negative xi", {
  params <- c(-.5, 1, .5)
  parametrisation <- "standard"
  target <- "transform"
  translation <- ParametrisationTranslator(
    params = params,
    parametrisation = parametrisation,
    target = target)
  expect_equal(translation, c(1., 1. + params[3], params[3]))
})


test_that("standard to noven - positive xi", {
  params <- c(.5, 1, .5)
  parametrisation <- "standard"
  target <- "noven"
  translation <- ParametrisationTranslator(
    params = params,
    parametrisation = parametrisation,
    target = target)
  expect_equal(translation, c(1./params[1], params[2]/params[1] - params[3], params[3]))
})

test_that("standard to noven - negative xi", {
  params <- c(-.5, 1, .5)
  parametrisation <- "standard"
  target <- "noven"
  translation <- ParametrisationTranslator(
    params = params,
    parametrisation = parametrisation,
    target = target)
  expect_equal(translation, c(1./params[1], params[2]/abs(params[1]) - params[3], params[3]))
})

test_that("noven to standard - positive xi", {
  params <- c(15., 5., .5)
  parametrisation <- "noven"
  target <- "standard"
  translation <- ParametrisationTranslator(
    params = params,
    parametrisation = parametrisation,
    target = target)
  expect_equal(translation, c(1./params[1], (params[2] + params[3]) / params[1], params[3]))
})

test_that("noven to standard - negative xi", {
  params <-  c(-15., 5., .5)
  parametrisation <- "noven"
  target <- "standard"
  translation <- ParametrisationTranslator(
    params = params,
    parametrisation = parametrisation,
    target = target)
  expect_equal(translation, c(1./params[1], (params[2] + params[3]) / abs(params[1]), params[3]))
})


test_that("wrong parametrisation", {
  params <-  c(-15., 5., .5)
  parametrisation <- "transform"
  target <- "standard"
  expect_error(ParametrisationTranslator(
    params = params,
    parametrisation = parametrisation,
    target = target))
})

