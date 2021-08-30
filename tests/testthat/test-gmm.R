
test_that("trawl objective", {
  n <- 1500
  print(print(list.files(path = ".")))
  print(print(list.files(path = "./../")))
  print(list.files(path = "./../../data/"))
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  depth <- 6
  data <- pollution_data[seq_len(n), test_column]
  pars_gpd <- evir::gpd(data, threshold = 0)$par.ses
  p_plus <- mean(data > 0)
  kappa <- 1 / p_plus - 1.
  pars <- c(pars_gpd, kappa)
  trawl_obj <- trawl_gmm$trawl_objective(data = data, depth = depth)
  trawl_obj_as_trawl_params <- trawl_obj(pars)

  trawl_values <- vapply(
    seq_len(7) / 20, function(x) trawl_obj_as_trawl_params(x), .1
  )
  testthat::expect_equal(which.min(trawl_values), 4)
})


test_that("trawl objective - grad", {
  n <- 1500
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]
  depth <- 4
  pars_gpd <- evir::gpd(pollution_data[, test_column], threshold = 0)$par.ses
  p_plus <- mean(pollution_data[, test_column] > 0)
  kappa <- 1 / p_plus - 1.
  pars <- c(pars_gpd, kappa)
  trawl_obj <- trawl_gmm$trawl_objective(
    data = pollution_data[, test_column], depth = depth
  )
  trawl_obj_as_trawl_params <- trawl_obj(pars)

  rho_vals <- seq_len(5) / 20 * .3
  trawl_grad_values <- lapply(
    X = rho_vals,
    FUN = function(x) pracma::grad(trawl_obj_as_trawl_params, x0 = x)
  )
  trawl_grad_values <- unlist(trawl_grad_values)

  testthat::expect_true(is.vector(trawl_grad_values))
  testthat::expect_equal(length(trawl_grad_values), length(rho_vals))
  testthat::expect_equal(which.min(abs(trawl_grad_values)), 4)
})

test_that("GMM objective - positive", {
  n <- 1500
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]
  init_guess_bds <- get_initial_guess_and_bounds(data = data)

  for (dp in c(3, 8)) {
    gmm_obj <- trawl_gmm$full_gmm_objective(data = data, depth = dp)
    gmm_vals <- gmm_obj(c(init_guess_bds$init_guess, .2))
    testthat::expect_true(gmm_vals > 0.0)
    testthat::expect_equal(length(gmm_vals), 1)
    testthat::expect_true(gmm_vals >= 0.0)
  }
})

test_that("Two-stage GMM objective - time & value", {
  n <- 1500
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]
  init_guess_bds <- get_initial_guess_and_bounds(data = data)
  two_step_gmm_obj <- trawl_gmm$two_stage_gmm_objective(data = data, depth = 3)

  start_time <- Sys.time()
  trawl_param_value <- optim(
    par = c(.2), fn = two_step_gmm_obj,
    method = "L-BFGS-B", lower = c(.001), upper = c(2)
  )$par
  time_delta <- Sys.time() - start_time

  testthat::expect_equal(trawl_param_value, .15, tolerance = 5e-2)
  testthat::expect_true(time_delta < 60)
})

test_that("Two-stage GMM objective - score", {
  n <- 3000
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]
  init_guess_bds <- get_initial_guess_and_bounds(data = data)
  max_length <- 50
  depth <- 3

  full_gmm_score <- trawl_gmm$trawl_gmm_score(
    params = c(init_guess_bds$init_guess, .15),
    depth = depth, type = "exp", max_length = max_length
  )
  score <- full_gmm_score(data)
  testthat::expect_false(any(vapply(score, function(score_per_depth) {
    any(is.na(score_per_depth))
  }, T)))
  testthat::expect_false(
    any(vapply(
      score, function(score_per_depth) any(is.infinite(score_per_depth)), T
    ))
  )
  dims_ground_vals <- vapply(
    seq_len(depth), function(i) c(max_length - i, 4), c(1, 1)
  )
  dim_score <- vapply(score, dim, c(1, 1))
  testthat::expect_equal(dim_score, dims_ground_vals)
})


test_that("Two-stage GMM objective - HAC full", {
  n <- 1500
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]
  init_guess_bds <- get_initial_guess_and_bounds(data = data)

  max_length <- 50
  depth <- 4
  k_max <- 5

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 10
  )
  hac_full <- trawl_gmm$trawl_gmm_hac(
    data = data, params = c(init_guess_bds$init_guess, i_guess),
    depth = depth, type = "exp",
    max_length = max_length, k = k_max
  )
  testthat::expect_true(Matrix::det(hac_full) > 0)
  testthat::expect_false(any(is.na(hac_full)))
  testthat::expect_false(any(is.infinite(hac_full)))
})

test_that("Two-stage GMM objective - HAC partial", {
  n <- 1500
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  data <- pollution_data[seq_len(n), test_column]
  init_guess_bds <- get_initial_guess_and_bounds(data = data)
  max_length <- 150
  depth <- 4
  k_max <- 5

  hac_partial <- trawl_gmm$trawl_gmm_hac_partial(
    data = data, params = c(init_guess_bds$init_guess, .15),
    depth = depth, type = "exp", max_length = max_length, k = k_max
  )
  testthat::expect_true(hac_partial > 0)
  testthat::expect_false(is.na(hac_partial))
  testthat::expect_false(is.infinite(hac_partial))
})

test_that("Two-stage - Variance", {
  n <- 1500
  pollution_data <- read.csv("./../../data/short_pollution_data.csv.gz")
  test_column <- 2
  depth <- 4
  max_depth <- 25
  data <- pollution_data[seq_len(n), test_column]

  i_guess <- pairwise_likelihood$init_guess(
    data = data, depth = depth, n_trials = 20
  )
  i_guess_model <- composite_marginal_mle(data)
  # init vals are c(-0.009792636, 0.3141497, 19.96388, 0.220771)
  ts_var <- trawl_gmm$two_stage_variance(
    data = data, params = c(i_guess_model, i_guess),
    depth = depth, type = "exp", max_length = max_depth
  )
  testthat::expect_true(ts_var > 0)
  testthat::expect_false(is.na(ts_var > 0))
})
