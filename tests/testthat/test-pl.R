
test_that("pair_pdf_constructor", {
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )
  testthat::expect_equal(
    pdf_constructor(xs = c(0, 1), h = 1), pdf_constructor(xs = c(1, 0), h = 1)
  )
  testthat::expect_equal(
    pdf_constructor(xs = c(0, 0), h = 1), pdf_constructor(xs = c(0, 0), h = 1)
  )
  testthat::expect_equal(
    pdf_constructor(xs = c(1, 1), h = 1), pdf_constructor(xs = c(1, 1), h = 1)
  )
})

test_that("pl_constructor - not parallel", {
  n <- 1000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- n
  params <- c(.1, 1., 19, .2)
  depth <- 3
  data <- pollution_data[seq_len(max_length), test_column]

  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )

  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params, depth = depth, pair_likehood = pdf_constructor
  )

  testthat::expect_false(is.na(pl_constructor(data)))
})

test_that("pl_constructor - parallel", {
  n <- 1000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- n
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )
  data <- pollution_data[seq_len(max_length), test_column]

  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(pmax(cores - 1, 1))
  parallel::clusterExport(
    cl, c(
      "transformation_map_inverse", "transformation_map",
      "transformation_jacobian", "parametrisation_translator",
      "pairwise_likelihood", get_trawl_envs_list()
    )
  )

  depth <- 3
  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params, depth = depth, pair_likehood = pdf_constructor, cl = cl
  )
  pl_lik <- pl_constructor(data)
  testthat::expect_equal(length(pl_lik), 1)
  testthat::expect_false(is.null(pl_lik))
  testthat::expect_false(any(is.na(pl_lik)))
  parallel::stopCluster(cl) # release resources
})

test_that("pl_constructor - parallel vs not parallel", {
  time_divisor <- 1e6
  n <- 3000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- n
  params <- c(.1, 1., 19, .2)
  depth <- 15
  data <- pollution_data[seq_len(max_length), test_column]

  # test times
  test_repeat <- 5

  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )

  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params, depth = depth, pair_likehood = pdf_constructor
  )

  res_no_parallel <- pl_constructor(data)
  no_parallel_times <- microbenchmark::microbenchmark(
    function() pl_constructor(data),
    times = test_repeat
  )$time / time_divisor
  no_parallel_times <- mean(no_parallel_times)

  cores <- parallel::detectCores()
  cl <- parallel::makeCluster(max(cores - 1, 1))
  parallel::clusterExport(
    cl, c(
      "transformation_map_inverse", "transformation_map",
      "transformation_jacobian", "parametrisation_translator",
      "pairwise_likelihood", get_trawl_envs_list()
    )
  )

  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params, depth = depth, pair_likehood = pdf_constructor, cl = cl
  )
  res_parallel <- pl_constructor(data)
  parallel_times <- microbenchmark::microbenchmark(
    function() pl_constructor(data),
    times = test_repeat
  )$time / time_divisor
  parallel_times <- mean(parallel_times)
  testthat::expect_equal(res_parallel, res_no_parallel, tolerance = 1e-3)
  testthat::expect_lte(parallel_times / no_parallel_times, 1.5)
  parallel::stopCluster(cl) # release resources
})

test_that("pl_constructor - PL initial guess", {
  time_divisor <- 1e6
  n <- 3000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- n
  depth <- 5
  data <- pollution_data[seq_len(max_length), test_column]

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 10
  )
  testthat::expect_equal(i_guess, .10, tolerance = 2e-2)
})

test_that("pl_constructor - PL as function of rho - convex", {
  time_divisor <- 1e6

  n <- 3000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- n
  depth <- 5

  cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(cores - 1)
  parallel::clusterExport(
    cl, c(
      "transformation_map_inverse", "transformation_map",
      "transformation_jacobian", "parametrisation_translator",
      "pairwise_likelihood", get_trawl_envs_list()
    )
  )

  pl_fn <- pairwise_likelihood$two_stage_trawl_pl(
    data = pollution_data[seq_len(max_length), test_column],
    depth = depth, cl = cl
  )
  rho_vals <- seq_len(10) / 30
  pl_rho_values <- vapply(rho_vals, pl_fn, 1.0)
  which_positive <- which(diff(pl_rho_values) > 0)
  which_negative <- which(diff(pl_rho_values) < 0)

  testthat::expect_true(length(which_positive) > 0)
  testthat::expect_true(length(which_positive) > 0)
  testthat::expect_true(max(which_negative) < min(which_positive)) # convexity
})

test_that("trawl_pl_hac", {
  time_divisor <- 1e6

  n <- 3000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- 1000
  depth <- 5

  data <- pollution_data[seq_len(max_length), test_column]
  k_max <- 20
  max_length <- 500

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 10
  )
  i_guess_model <- get_initial_guess_and_bounds(
    data = data, max_length = max_length
  )
  hac_full <- pairwise_likelihood$trawl_pl_hac(
    data = data, params = c(i_guess_model$init_guess, i_guess),
    depth = depth, k = k_max, type = "exp", max_length = max_length
  )

  testthat::expect_true(Matrix::det(hac_full) > 0)
  testthat::expect_false(any(is.na(hac_full)))
  testthat::expect_false(any(is.infinite(hac_full)))
})

test_that("trawl_pl_hac_partial", {
  time_divisor <- 1e6

  n <- 3000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- 1000
  depth <- 5

  data <- pollution_data[seq_len(max_length), test_column]

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 10
  )
  i_guess_model <- get_initial_guess_and_bounds(
    data = data, max_length = max_length
  )
  hac_partial <- pairwise_likelihood$trawl_pl_hac_partial(
    data = data, params = c(i_guess_model$init_guess, i_guess),
    depth = depth, k = 8, type = "exp", max_length = 500
  )
  testthat::expect_true(hac_partial > 0)
  testthat::expect_false(is.na(hac_partial))
  testthat::expect_false(is.infinite(hac_partial))
})


test_that("PLHessian", {
  time_divisor <- 1e6
  n <- 3000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- 1000
  depth <- 4

  data <- pollution_data[seq_len(max_length), test_column]

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 10
  )
  i_guess_model <- get_initial_guess_and_bounds(
    data = data, max_length = max_length
  )
  pl_hessian <- pairwise_likelihood$trawl_pl_hessian(
    params = c(i_guess_model$init_guess, i_guess),
    depth = depth, type = "exp", max_length = 200
  )
  pl_hess <- pl_hessian(data[1:50])
  lapply(pl_hess, function(x) {
    testthat::expect_false(any(is.na(x)))
    testthat::expect_false(any(is.infinite(x)))
  })

  ts_var <- pairwise_likelihood$two_stage_variance(
    data = data, params = c(i_guess_model$init_guess, i_guess),
    depth = depth, type = "exp", max_length = 100
  )
  testthat::expect_true(ts_var > 0)
  testthat::expect_false(is.na(ts_var))
  testthat::expect_false(is.infinite(ts_var))
})


test_that("Two-stage - Variance", {
  n <- 3000
  pollution_data <- read.csv("../../data/short_pollution_data.csv.gz")
  test_column <- 2
  max_length <- 1000
  depth <- 4

  data <- pollution_data[seq_len(max_length), test_column]

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 20
  )
  i_guess_model <- composite_marginal_mle(data)
  # init should be c(-0.009792636, 0.3141497, 19.96388, 0.220771)
  ts_var <- pairwise_likelihood$two_stage_variance(
    data = data, params = c(i_guess_model, i_guess),
    depth = depth, type = "exp", max_length = 200
  )
  testthat::expect_true(ts_var > 0)
  testthat::expect_false(is.na(ts_var))
  testthat::expect_false(is.infinite(ts_var))
})
