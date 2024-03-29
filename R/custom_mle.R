
# Generates a GPD initial guess and parameter bounds for "L-BFGS-B".
get_initial_guess_and_bounds <- function(data, max_length = 20000,
                                         minus_mult = 0.5,
                                         plus_multiplier = 1.5) {
  length_data <- min(max_length, length(data))
  data <- data[1:length_data]

  invisible(utils::capture.output(
    init_guess_vals <- as.numeric(
      evir::gpd(data, threshold = 0, method = "pwm")$par.ests
    )
  )) # (xi, sigma)

  p_above_zero <- length(which(data > 0)) / length(data)
  kap <- (1 - p_above_zero) / p_above_zero
  init_guess_vals <- c(init_guess_vals, kap)
  minus_init <- minus_mult * init_guess_vals
  plus_init <- plus_multiplier * init_guess_vals
  bds <- cbind(minus_init, plus_init)
  lower <- apply(bds, 1, min)
  upper <- apply(bds, 1, max)

  return(list(init_guess = init_guess_vals, lower = lower, upper = upper))
}


#' Likelihood function of parameters for fixed data.
#' @param data Data on which to compute likelihood.
#' @return Custom likelihood function of parameters.
#' @examples
#' set.seed(42)
#' xi <- 1.
#' kappa <- 9.
#' sigma <- 1.
#' n <- 1000
#' p_zero <- 1 - 1 / (1 + kappa)
#' zeroes <- stats::runif(n = n, min = 0, max = 1)
#' test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
#' test_samples[which(zeroes < p_zero)] <- 0.0
#' cm_mle <- custom_likelihood(data = test_samples)
#' @export
#' @importFrom evir dgpd
custom_likelihood <- function(data) {
  assertthat::assert_that(is.vector(data))

  return(function(par) {
    data_for_mle <- data
    kap <- par[3]
    p_zero <- 1 - 1. / (1. + kap)

    invisible(utils::capture.output(
      like <- evir::dgpd(
        x = data_for_mle[data_for_mle > 0.0],
        beta = par[2], xi = par[1], mu = 0.0
      ),
      type = "message"
    ))
    like <- like * (1 - p_zero)
    like <- like[!is.na(like)] + 1e-7
    like <- sum(log(like))

    log.like <- like + length(which(data_for_mle == 0.0)) * log(p_zero)
    log.like <- min(log.like, 1e9)
    log.like <- max(log.like, -1e9)

    return(-log.like)
  })
}


custom_marginal_mle <- function(data, parametrisation = "standard") {
  get_init_guess <- get_initial_guess_and_bounds(data)
  init_guess <- get_init_guess$init_guess
  lower <- get_init_guess$lower
  upper <- get_init_guess$upper
  fn_mle <- custom_likelihood(data = data)

  res <- stats::optim(
    par = init_guess[1:3], fn_mle, method = "L-BFGS-B",
    lower = lower, upper = upper
  )$par
  return(res)
}


composite_likelihood <- function(data) {
  # implements custom likelihood with jacobian
  # returns a function of parameters
  return(function(par) {
    data_for_mle <- data

    trf_data <- transformation_map_inverse(data_for_mle, params_std = par)

    # needs to inject transformed data here
    custom_lk <- custom_likelihood(data = trf_data)
    # computes GPD(1, 1+kappa)
    trf_lk <- custom_lk(par = c(1., 1. + par[3], par[3]))

    jacob <- transformation_jacobian(params_std = par)
    jacob_vals <- jacob(data_for_mle)
    jacob_vals <- -sum(log(jacob_vals[!is.na(jacob_vals)] + 1e-7))

    composite_lk <- jacob_vals + trf_lk
    return(composite_lk)
  })
}


composite_marginal_mle <- function(data) {
  get_init_guess <- get_initial_guess_and_bounds(data)
  init_guess <- get_init_guess$init_guess
  lower <- get_init_guess$lower
  upper <- get_init_guess$upper
  fn_mle <- composite_likelihood(data = data)

  res <- stats::optim(
    par = init_guess[1:3], fn_mle, method = "L-BFGS-B",
    lower = lower, upper = upper, control = list()
  )$par
  return(res)
}


composite_likelihood_score <- function(params, max_length = 100) {
  return(
    function(data) {
      n_row <- min(max_length, length(data))
      data <- data[1:n_row]
      grad_composite <- t(vapply(data, function(x) {
        composite_lk <- composite_likelihood(x)
        return(pracma::grad(composite_lk, x0 = params))
      }, rep(0, length(params))))
      return(grad_composite)
    }
  ) # data_length x length(params)
}

composite_likelihood_hessian <- function(params, max_length = 100) {
  return(
    function(data) {
      n_row <- min(max_length, length(data))
      data <- data[1:n_row]
      hess_composite <- t(vapply(data, function(x) {
        composite_lk <- composite_likelihood(x)
        return(pracma::hessian(composite_lk, x0 = params))
      }, rep(0, length(params)^2)))
      return(hess_composite)
    }
  ) # data_length x length(params)^2
}

#' HAC estimator of autocovariance for the composite marginal.
#' @param data Data on which to compute the autocovariance.
#' @param params Composite likelihood params.
#' @param k Maximum lag.
#' @param max_length Maximum length of data used to compute the score functions.
#' @param near_pd Boolean, nearest positive definite matrix from HAC estimate.
#' @return Composite marginal HAC autocovariance.
#' @examples
#' xi <- 1.
#' kappa <- 9.
#' sigma <- 1.
#' n <- 100
#' p_zero <- 1 - 1. / (1. + kappa)
#' test_column <- 2
#'
#' set.seed(42)
#' zeroes <- stats::runif(n = n, min = 0, max = 1)
#' test_samples <- evir::rgpd(n = n, xi = xi, mu = 0, beta = sigma)
#' test_samples[which(zeroes < p_zero)] <- 0.0
#' composite_marginal_hac(
#'   data = test_samples, params = c(xi, sigma, kappa),
#'   k = 3, max_length = 50
#' )
#' @export
composite_marginal_hac <- function(data, params, k = 10,
                                   max_length = 100, near_pd = F) {
  lk_score <- composite_likelihood_score(params, max_length)
  score_acf <- autocovariance_matrix(lk_score(data), k)
  return(make_hac(score_acf, near_pd))
}
