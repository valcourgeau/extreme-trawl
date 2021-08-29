

test_that("trawl_simulation - NAs, Gamma and shape", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 1500
  vd <- 100
  type <- "exp"
  trawl_parameter <- 0.1

  trawl_sims <- trawl_simulation(
    alpha = 1., beta = 1., trawl_parameter = trawl_parameter, n = n,
    vanishing_depth = vd, type = type, parallel = F
  )
  testthat::expect_equal(length(trawl_sims), n)
  testthat::expect_false(any(is.na(trawl_sims)))
  invisible(capture.output(fit_gam <- fitdistrplus::fitdist(
    data = trawl_sims, distr = "gamma", method = "mle"
  )))
  testthat::expect_equal(unname(fit_gam$estimate), c(1, 1), tolerance = 0.05)
})

test_that("exceedances_simulation - simple", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 3000
  vd <- 100
  type <- "exp"

  exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd, type = type
  )
  # returns (exceedances, latent)
  gpd_fit_sim <- evir::gpd(exc$exceedances, threshold = 0.0, method = "ml")
  par_vals <- unname(gpd_fit_sim$par.ests)
  tol_vals <- unname(gpd_fit_sim$par.ses)
  testthat::expect_equal(par_vals[1], params[1], tolerance = 1.96 * tol_vals[1])
  testthat::expect_equal(par_vals[2], params[2], tolerance = 1.96 * tol_vals[2])
  testthat::expect_equal(mean(exc$exceedances > 0), .05, tolerance = .05)

  acf_vals <- acf(exc$exceedances, lag.max = 10, plot = F)$acf[, , 1]
  truth <- trawl_autocorrelation$acf_trawl_collection(
    h = c(.01, 1:10), alpha = 1, beta = 1, kappa = params[3], rho = params[4],
    cov = F, type = type, delta = .1, end_seq = 100
  )
})

test_that("exceedances_simulation - cross", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 2500
  vd <- 50
  type <- "exp"

  testthat::expect_warning(exceedances_simulation(
    params = params, n = n, vanishing_depth = vd,
    type = type, m = 2, algo = "cross"
  ))
  suppressWarnings(exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd,
    type = type, m = 50, algo = "cross"
  ))
  testthat::expect_false(any(is.na(exc$exceedances)))
  testthat::expect_true(any(is.na(exc$latent)))

  testthat::expect_equal(
    mean(exc$exceedances > 0, na.rm = T), .05, tolerance = .07
  )
})

test_that("exceedances_simulation - corr unif", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 2000
  vd <- 10
  type <- "exp"

  truth <- trawl_autocorrelation$acf_trawl_collection(
    h = c(.01, 1:10), alpha = 1, beta = 1, kappa = params[3], rho = params[4],
    cov = F, type = type, delta = .1, end_seq = 100
  )

  vd_list <- c(5, 10, 20, 50, 70, 100, 200)
  vd_error <- rep(0, length(vd_list))
  cove_error <- rep(0, length(vd_list))
  i <- 1
  for (vd in vd_list) {
    cove <- print_vanishing_coverage(
      trawl_parameter = params[4], vanishing_depth = vd,
      type = type, get_value = T
    )
    exc <- exceedances_simulation(
      params = params, n = n, vanishing_depth = vd, type = type,
      m = 3, algo = "corr_unif"
    )
    acf_vals <- acf(exc$exceedances, plot = F, lag.max = 10)$acf[, , 1]
    vd_error[i] <- sqrt(sum((acf_vals - truth)^2))
    cove_error[i] <- cove
    i <- i + 1
    testthat::expect_true(sum(diff(acf_vals) < 0) > 5)
  }
  # WE CAN PLOT VD LIST AGAINST COV ERROR & VD ERROR
})

test_that("exceedances_simulation - dynamic latent", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 3000
  vd <- 100
  type <- "exp"

  exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd, type = type,
    m = 3, algo = "dynamic_latent"
  )
  acf_vals <- acf(exc$exceedances, plot = F, lag.max = 10)$acf[, , 1]
  testthat::expect_true(sum(diff(acf_vals) < 0) > 5)
  testthat::expect_equal(
    mean(exc$exceedances > 0, na.rm = T), .05,
    tolerance = .1
  )
})

test_that("exceedances_simulation - dynamic uniform", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 3000
  vd <- 100
  type <- "exp"

  exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd, type = type,
    m = 3, algo = "dynamic_uniform"
  )
  acf_vals <- acf(exc$exceedances, plot = F, lag.max = 10)$acf[, , 1]
  testthat::expect_true(sum(diff(acf_vals) < 0) > 5)
  testthat::expect_equal(
    mean(exc$exceedances > 0, na.rm = T), .05,
    tolerance = .1
  )
})
