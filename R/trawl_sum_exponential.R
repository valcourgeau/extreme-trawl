sum_exponential_trawl <- new.env()

sum_exponential_trawl$sum_exp <- function(param, h, trawl_f) {
  stopifnot(all(param > 0))
  n_sup <- length(param)
  exp_trawl_params <- param
  trawl_val <- vapply(exp_trawl_params, function(par) {
    trawl_f(param = par, h)
  }, rep(1.0, length(h)))

  trawl_val <- trawl_val * matrix(
    rep(param / sum(param), length(h)),
    ncol = n_sup, byrow = T
  )
  return(rowSums(trawl_val))
}

sum_exponential_trawl$trawl_b_one <- function(param, h) {
  return(
    sum_exponential_trawl$sum_exp(param, h, exponential_trawl$trawl_b_one)
  )
}

sum_exponential_trawl$trawl_b_two <- function(param, h) {
  return(
    sum_exponential_trawl$sum_exp(param, h, exponential_trawl$trawl_b_two)
  )
}

sum_exponential_trawl$trawl_b_three <- function(param, h) {
  return(
    sum_exponential_trawl$sum_exp(param, h, exponential_trawl$trawl_b_three)
  )
}

sum_exponential_trawl$config <- function() {
  n_params <- 4
  return(
    list(
      "n_params" = n_params,
      "lower" = rep(1e-3, n_params),
      "upper" = rep(1.0, n_params)
    )
  )
}
