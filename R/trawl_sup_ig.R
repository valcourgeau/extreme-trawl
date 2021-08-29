sup_ig_trawl <- new.env()

sup_ig_trawl$leb_a <- function(param) {
  # param is (gamma, delta)
  stopifnot(length(param) == 2)
  stopifnot(all(param > 0))
  return(param[1] / param[2])
}

sup_ig_trawl$trawl_b_one <- function(param, h) {
  stopifnot(length(param) == 2)
  stopifnot(all(param > 0))
  a_total <- sup_ig_trawl$leb_a(param)
  log_scaling_factor <- param[1] * param[2] * (1 - sqrt(1 + 2 * h / param[1]^2))
  return(a_total * (1 - exp(log_scaling_factor)))
}

sup_ig_trawl$trawl_b_two <- function(param, h) {
  stopifnot(length(param) == 2)
  stopifnot(all(param > 0))
  a_total <- sup_ig_trawl$leb_a(param)
  log_scaling_factor <- param[1] * param[2] * (1 - sqrt(1 + 2 * h / param[1]^2))
  return(a_total * exp(log_scaling_factor))
}

sup_ig_trawl$trawl_b_three <- function(param, h) {
  stopifnot(length(param) == 2)
  stopifnot(all(param > 0))
  return(sup_ig_trawl$trawl_b_one(param, h))
}

sup_ig_trawl$config <- function() {
  return(list(n_params = 2, lower = rep(1e-5, 2), upper = rep(10.0, 2)))
}
