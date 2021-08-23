test_that("trawl objective", {
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  pars_gpd <- evir::gpd(pollution_data[, test_column], threshold = 0)$par.ses
  p_plus <- mean(pollution_data[, test_column] > 0)
  kappa <- 1 / p_plus - 1.
  pars <- c(pars_gpd, kappa)
  trawl_obj <- trawl_gmm$trawl_objective(
    data = pollution_data[, test_column],
    depth = 10
  )
  trawl_obj_as_trawl_params <- trawl_obj(pars)

  trawl_values <- vapply(1:10 / 10, function(x) {
    trawl_obj_as_trawl_params(x)
  }, .1)
  testthat::expect_equal(which.min(trawl_values), 2)
})


test_that("trawl objective - grad", {
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  pars_gpd <- evir::gpd(pollution_data[, test_column], threshold = 0)$par.ses
  p_plus <- mean(pollution_data[, test_column] > 0)
  kappa <- 1 / p_plus - 1.
  pars <- c(pars_gpd, kappa)
  trawl_obj <- trawl_gmm$trawl_objective(
    data = pollution_data[, test_column],
    depth = 10
  )
  trawl_obj_as_trawl_params <- trawl_obj(pars)

  trawl_grad_values <- lapply(as.list(1:5 / 5 * .3), function(x) {
    pracma::grad(trawl_obj_as_trawl_params, x0 = x)
  })
  trawl_grad_values <- unlist(trawl_grad_values)
  testthat::expect_equal(which.min(abs(trawl_grad_values)), 3)
})

test_that("GMM objective - positive", {
  max_length <- 30000
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  init_guess_bds <- get_initial_guess_and_bounds(
    data = pollution_data[1:max_length, test_column]
  )

  gmm_obj <- trawl_gmm$full_gmm_objective(
    data = pollution_data[, test_column],
    depth = 10
  )

  testthat::expect_true(gmm_obj(c(init_guess_bds$init_guess, .2)) > 0.0)
})

test_that("GMM objective", {
  max_length <- 30000
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  init_guess_bds <- get_initial_guess_and_bounds(
    data = pollution_data[1:max_length, test_column]
  )

  gmm_obj <- trawl_gmm$full_gmm_objective(
    data = pollution_data[, test_column],
    depth = 3
  )

  print(microbenchmark::microbenchmark(
    gmm_obj(c(init_guess_bds$init_guess, .2)),
    times = 5
  ))
  print(optim(
    par = c(init_guess_bds$init_guess, .2),
    fn = gmm_obj, method = "L-BFGS-B",
    lower = c(init_guess_bds$lower, .001),
    upper = c(init_guess_bds$upper, 2),
    control = list(trace = 3)
  )$par)

  testthat::expect_equal(T, T)
})

test_that("Two-stage GMM objective - time & value", {
  max_length <- 30000
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  init_guess_bds <- get_initial_guess_and_bounds(
    data = pollution_data[1:max_length, test_column]
  )

  two_step_gmm_obj <- trawl_gmm$two_stage_gmm_objective(
    data = pollution_data[, test_column],
    depth = 10
  )

  start_time <- Sys.time()
  trawl_param_value <- optim(
    par = c(.2),
    fn = two_step_gmm_obj, method = "L-BFGS-B",
    lower = c(.001),
    upper = c(2),
  )$par
  time_delta <- Sys.time() - start_time

  testthat::expect_equal(trawl_param_value, .15, tolerance = 5e-2)
  testthat::expect_true(time_delta < 60)
})

test_that("Two-stage GMM objective - score", {
  max_length <- 30000
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  data <- pollution_data[1:max_length, test_column]
  init_guess_bds <- get_initial_guess_and_bounds(
    data = data
  )
  max_length <- 100
  depth <- 3

  full_gmm_score <- trawl_gmm$trawl_gmm_score(
    params = c(init_guess_bds$init_guess, .15),
    depth = depth, type = "exp", max_length = max_length
  )
  score <- full_gmm_score(data)
  testthat::expect_false(any(vapply(score, function(score_per_depth) {
    any(is.na(score_per_depth))
  }, T)))
  testthat::expect_false(any(vapply(score, function(score_per_depth) {
    any(is.infinite(score_per_depth))
  }, T)))
})

test_that("Two-stage GMM objective - HAC full", {
  max_length <- 30000
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  data <- pollution_data[1:max_length, test_column]
  init_guess_bds <- get_initial_guess_and_bounds(
    data = data
  )
  max_length <- 500
  depth <- 12

  k_max <- 20

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 10
  )
  hac_full <- trawl_gmm$trawl_gmm_hac(
    data = data,
    params = c(init_guess_bds$init_guess, i_guess),
    depth = depth, type = "exp",
    max_length = max_length, k = k_max
  )
  print(hac_full)
  testthat::expect_true(Matrix::det(hac_full) > 0)
  testthat::expect_false(any(vapply(hac_full, is.na, T)))
  testthat::expect_false(any(vapply(hac_full, is.infinite, T)))
})

test_that("Two-stage GMM objective - HAC partial", {
  max_length <- 30000
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  data <- pollution_data[1:max_length, test_column]
  init_guess_bds <- get_initial_guess_and_bounds(
    data = data
  )
  max_length <- 500
  depth <- 3

  k_max <- 10

  hac_partial <- trawl_gmm$trawl_gmm_hac_partial(
    data = data,
    params = c(init_guess_bds$init_guess, .15), depth = depth, type = "exp",
    max_length = max_length, k = k_max
  )
  testthat::expect_true(hac_partial > 0)
  testthat::expect_false(is.na(hac_partial))
  testthat::expect_false(is.infinite(hac_partial))
})

test_that("Two-stage - Variance", {
  pollution_data <- read.csv("data/clean_pollution_data.csv")
  test_column <- 2
  max_length <- 20000
  depth <- 4

  data <- pollution_data[1:max_length, test_column]

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 20
  )
  i_guess_model <- composite_marginal_mle(data)
  # init vals are c(-0.009792636, 0.3141497, 19.96388, 0.220771)
  ts_var <- trawl_gmm$two_stage_variance(
    data = data, params = c(i_guess_model, i_guess),
    depth = depth, type = "exp", max_length = 200
  )
  cat("ts_var", ts_var, "\n")
})
