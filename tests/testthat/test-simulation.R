test_that("exceedances_simulation - simple", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 10000
  vd <- 100
  type <- "exp"

  exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd, type = type
  )
  # returns (exceedances, latent)
  gpd_fit_sim <- evir::gpd(exc$exceedances, threshold = 0.0, method = "ml")

  testthat::expect_equal(
    unname(gpd_fit_sim$par.ests)[1],
    params[1],
    tolerance = 1.96 * unname(gpd_fit_sim$par.ses)[1]
  )
  testthat::expect_equal(
    unname(gpd_fit_sim$par.ests)[2],
    params[2],
    tolerance = 1.96 * unname(gpd_fit_sim$par.ses)[2]
  )

  testthat::expect_equal(mean(exc$exceedances > 0), .05, tolerance = 8e-3)
  acf_vals <- acf(exc$exceedances, lag.max = 10, plot = F)$acf[, , 1]
  truth <- trawl_autocorrelation$acf_trawl_collection(
    h = c(.01, 1:10), alpha = 1, beta = 1, kappa = params[3], rho = params[4],
    cov = F, type = type, delta = .1, end_seq = 100
  )
})

test_that("exceedances_simulation - cross", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 5000
  vd <- 50
  type <- "exp"

  exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd,
    type = type, m = 3, algo = "cross"
  )
  # returns (exceedances, latent)
  gpd_fit_sim <- evir::gpd(exc$exceedances, threshold = 0.0, method = "ml")

  testthat::expect_equal(
    unname(gpd_fit_sim$par.ests)[1],
    params[1],
    tolerance = 1.96 * unname(gpd_fit_sim$par.ses)[1]
  )
  testthat::expect_equal(
    unname(gpd_fit_sim$par.ests)[2],
    params[2],
    tolerance = 1.96 * unname(gpd_fit_sim$par.ses)[2]
  )

  testthat::expect_equal(mean(exc$exceedances > 0), .05, tolerance = 8e-3)
  print(as.vector(acf(exc$exceedances, lag.max = 10, plot = F)$acf))
  print(trawl_autocorrelation$acf_trawl_collection(
    h = c(.01, 1:10), alpha = 1, beta = 1, kappa = params[3], rho = params[4],
    cov = F, type = type, delta = .1, end_seq = 100
  ))
})

test_that("exceedances_simulation - corr unif", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 5000
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
    acf_vals <- acf(
      exc$exceedances[1:(n - 2 * vd)],
      plot = F, lag.max = 10
    )$acf[, , 1]
    vd_error[i] <- sum((acf_vals - truth)^2)
    cove_error[i] <- cove
    i <- i + 1
  }

  plot(vd_list, vd_error, main = "Vanishing Depth Error", ylim = c(0, 3))
  lines(vd_list, cove_error, col = "red")
})

test_that("exceedances_simulation - dynamic latent", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 5000
  vd <- 10
  type <- "exp"

  truth <- trawl_autocorrelation$acf_trawl_collection(
    h = c(.01, 1:10), alpha = 1, beta = 1, kappa = params[3], rho = params[4],
    cov = F, type = type, delta = .1, end_seq = 100
  )

  vd <- 100
  exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd, type = type,
    m = 3, algo = "dynamic_latent"
  )
  cat(
    "Prob of positive",
    sum(exc$exceedances > 0) / length(exc$exceedances), "\n"
  )
  acf(exc$exceedances, main = "dynamic latent acf")
  lines(0:10, truth, col = "red")
})

test_that("exceedances_simulation - dynamic uniform", {
  set.seed(42)
  params <- c(.1, 1., 19, .05)
  n <- 5000
  vd <- 10
  type <- "exp"

  truth <- trawl_autocorrelation$acf_trawl_collection(
    h = c(.01, 1:10), alpha = 1, beta = 1, kappa = params[3], rho = params[4],
    cov = F, type = type, delta = .1, end_seq = 100
  )

  vd <- 100
  exc <- exceedances_simulation(
    params = params, n = n, vanishing_depth = vd, type = type,
    m = 3, algo = "dynamic_uniform"
  )
  cat(
    "Prob of positive", sum(exc$exceedances > 0) / length(exc$exceedances), "\n"
  )
  acf(exc$exceedances, main = "dynamic uniform acf")
  lines(0:10, truth, col = "red")
})
