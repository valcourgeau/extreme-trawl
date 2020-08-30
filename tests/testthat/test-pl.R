test_that("PairPDFConstructor", {
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- PairwiseLikelihood$PairPDFConstructor(
    params=params, type='exp'
  )
  testthat::expect_equal(pdf_constructor(xs=c(0,1), h=1), pdf_constructor(xs=c(1,0), h=1))
})


test_that("PLConstructor - not parallel", {
  pollution_data <- read.csv('../../data/clean_pollution_data.csv')
  test_column <- 2
  max_length <- 1000
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- PairwiseLikelihood$PairPDFConstructor(
    params=params, type='exp'
  )

  depth <- 3
  pl_constructor <- PairwiseLikelihood$PLConstructor(
    params=params,
    depth=depth,
    pair_likehood=pdf_constructor
  )


  testthat::expect_false(is.na(pl_constructor(pollution_data[1:max_length, test_column])))
})

test_that("PLConstructor - parallel", {
  pollution_data <- read.csv('../../data/clean_pollution_data.csv')
  test_column <- 2
  max_length <- 1000
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- PairwiseLikelihood$PairPDFConstructor(
    params=params, type='exp'
  )

  cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(cores-1)
  parallel::clusterExport(
    cl, c('TransformationMapInverse',
          'TransformationMap',
          'TransformationJacobian',
          'ParametrisationTranslator',
          'PairwiseLikelihood',
          GetTrawlEnvsList()))

  depth <- 3
  pl_constructor <- PairwiseLikelihood$PLConstructor(
    params=params,
    depth=depth,
    pair_likehood=pdf_constructor,
    cl=cl
  )
  # parallel::stopCluster(cl)

  testthat::expect_false(is.na(pl_constructor(pollution_data[1:max_length, test_column])))
})

test_that("PLConstructor - parallel vs not parallel", {
  TIME_DIVISOR <- 1e6

  pollution_data <- read.csv('../../data/clean_pollution_data.csv')
  test_column <- 2
  max_length <- 20000
  params <- c(.1, 1., 19, .2)
  depth <- 3

  pdf_constructor <- PairwiseLikelihood$PairPDFConstructor(
    params=params, type='exp'
  )

  pl_constructor <- PairwiseLikelihood$PLConstructor(
    params=params,
    depth=depth,
    pair_likehood=pdf_constructor
  )

  res_no_parallel <- pl_constructor(pollution_data[1:max_length, test_column])
  no_parallel_times <- microbenchmark::microbenchmark(pl_constructor(pollution_data[1:max_length, test_column]), times=5)$time / TIME_DIVISOR

  cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(cores-1)
  parallel::clusterExport(
    cl, c('TransformationMapInverse',
          'TransformationMap',
          'TransformationJacobian',
          'ParametrisationTranslator',
          'PairwiseLikelihood',
          GetTrawlEnvsList()))

  pl_constructor <- PairwiseLikelihood$PLConstructor(
    params=params,
    depth=depth,
    pair_likehood=pdf_constructor,
    cl=cl
  )
  res_parallel <- pl_constructor(pollution_data[1:max_length, test_column])
  parallel_times <- microbenchmark::microbenchmark(pl_constructor(pollution_data[1:max_length, test_column]), times=5)$time / TIME_DIVISOR
  parallel::stopCluster(cl)

  cat('PL Parallel improvement: ', round(mean(no_parallel_times / parallel_times)*100, 0), '%\n')
  testthat::expect_equal(res_parallel, res_no_parallel, tolerance=1e-3)
  testthat::expect_equal(parallel_times / no_parallel_times, rep(0, length(no_parallel_times)), tolerance=.5)
})




