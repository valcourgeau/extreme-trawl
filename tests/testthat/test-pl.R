test_that("pair_pdf_constructor", {
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )
  testthat::expect_equal(
    pdf_constructor(xs = c(0, 1), h = 1),
    pdf_constructor(xs = c(1, 0), h = 1)
  )
})

test_that("pl_constructor - not parallel", {
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 1000
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )

  depth <- 3
  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params,
    depth = depth,
    pair_likehood = pdf_constructor
  )

  testthat::expect_false(
    is.na(pl_constructor(pollution_data[1:max_length, test_column]))
  )
})

test_that("pl_constructor - parallel", {
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 1000
  params <- c(.1, 1., 19, .2)
  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )

  cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(cores - 1)
  parallel::clusterExport(
    cl, c(
      "transformation_map_inverse",
      "transformation_map",
      "transformation_jacobian",
      "parametrisation_translator",
      "pairwise_likelihood",
      get_trawl_envs_list()
    )
  )

  depth <- 3
  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params,
    depth = depth,
    pair_likehood = pdf_constructor,
    cl = cl
  )
  parallel::stopCluster(cl)

  testthat::expect_false(
    is.na(pl_constructor(pollution_data[1:max_length, test_column]))
  )
})

test_that("pl_constructor - parallel vs not parallel", {
  time_divisor <- 1e6

  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  params <- c(.1, 1., 19, .2)
  depth <- 3

  pdf_constructor <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = "exp"
  )

  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params,
    depth = depth,
    pair_likehood = pdf_constructor
  )

  res_no_parallel <- pl_constructor(pollution_data[1:max_length, test_column])
  no_parallel_times <- microbenchmark::microbenchmark(
    pl_constructor(pollution_data[1:max_length, test_column]),
    times = 5
  )$time / time_divisor

  cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(cores - 1)
  parallel::clusterExport(
    cl, c(
      "transformation_map_inverse",
      "transformation_map",
      "transformation_jacobian",
      "parametrisation_translator",
      "pairwise_likelihood",
      get_trawl_envs_list()
    )
  )

  pl_constructor <- pairwise_likelihood$pl_constructor(
    params = params,
    depth = depth,
    pair_likehood = pdf_constructor,
    cl = cl
  )
  res_parallel <- pl_constructor(pollution_data[1:max_length, test_column])
  parallel_times <- microbenchmark::microbenchmark(
    pl_constructor(pollution_data[1:max_length, test_column]),
    times = 5
  )$time / time_divisor
  parallel::stopCluster(cl)

  cat(
    "\nPL Parallel improvement:",
    round(mean(no_parallel_times / parallel_times) * 100, 0),
    "%\n"
  )
  testthat::expect_equal(
    res_parallel, res_no_parallel,
    tolerance = 1e-3
  )
  testthat::expect_equal(
    parallel_times / no_parallel_times,
    rep(0, length(no_parallel_times)),
    tolerance = .5
  )
})


test_that("pl_constructor - PL initial guess", {
  time_divisor <- 1e6

  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  depth <- 10

  i_guess <- pairwise_likelihood$init_guess(
    data = pollution_data[1:max_length, test_column],
    depth = depth, n_trials = 10
  )
  testthat::expect_equal(i_guess, .15, tolerance = 1e-2)
})

test_that("pl_constructor - PL as function of rho - convex", {
  time_divisor <- 1e6

  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  depth <- 10

  cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(cores - 1)
  parallel::clusterExport(
    cl, c(
      "transformation_map_inverse",
      "transformation_map",
      "transformation_jacobian",
      "parametrisation_translator",
      "pairwise_likelihood",
      get_trawl_envs_list()
    )
  )

  pl_fn <- pairwise_likelihood$two_stage_trawl_pl(
    data = pollution_data[1:max_length, test_column], depth = depth, cl = cl
  )
  rho_vals <- 1:10 / 30
  pl_rho_values <- vapply(rho_vals, pl_fn, 1.0)
  which_positive <- which(diff(pl_rho_values) > 0)
  which_negative <- which(diff(pl_rho_values) < 0)

  plot(rho_vals, pl_rho_values, main = "pl as function of rho")

  testthat::expect_true(length(which_positive) > 0)
  testthat::expect_true(length(which_positive) > 0)
  testthat::expect_true(max(which_negative) < min(which_positive))
})


test_that("trawl_pl_hac", {
  time_divisor <- 1e6

  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  depth <- 12

  data <- pollution_data[1:max_length, test_column]
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
  print("hac_full")
  print(hac_full)
  testthat::expect_true(Matrix::det(hac_full) > 0)
  testthat::expect_false(any(vapply(hac_full, is.na, T)))
  testthat::expect_false(any(vapply(hac_full, is.infinite, T)))
})

test_that("trawl_pl_hac_partial", {
  time_divisor <- 1e6

  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  depth <- 12

  data <- pollution_data[1:max_length, test_column]

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
  cat("hac_partial", hac_partial, "\n")
  testthat::expect_true(hac_partial > 0)
  testthat::expect_false(is.na(hac_partial))
  testthat::expect_false(is.infinite(hac_partial))
})


test_that("PLHessian", {
  time_divisor <- 1e6

  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  depth <- 12

  data <- pollution_data[1:max_length, test_column]

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
  print("pl_hessian")
  pl_hess <- pl_hessian(data[1:50])
  print(matrix(pl_hess[[1]][3, ], 4, 4))

  ts_var <- pairwise_likelihood$two_stage_variance(
    data = data, params = c(i_guess_model$init_guess, i_guess),
    depth = depth, k = 10, type = "exp", max_length = 100
  )
})


test_that("Two-stage - Variance", {
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  depth <- 8

  data <- pollution_data[1:max_length, test_column]

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 20
  )
  i_guess_model <- composite_marginal_mle(data)
  # init should be c(-0.009792636, 0.3141497, 19.96388, 0.220771)
  ts_var <- pairwise_likelihood$two_stage_variance(
    data = data, params = c(i_guess_model, i_guess),
    depth = depth, type = "exp", max_length = 200
  )
  cat("ts_var", ts_var, "\n")
})
