gamma_trawl <- new.env()

gamma_trawl$leb_a <- function(param) {
  # param is (alpha, H)
  stopifnot(length(param) == 2)
  stopifnot(all(param > 0))
  return(param[1] / (param[2] - 1))
}

gamma_trawl$trawl_b_one <- function(param, h) {
  stopifnot(length(param) == 2)
  stopifnot(all(param > 0))
  a_total <- gamma_trawl$leb_a(param)
  return(a_total - gamma_trawl$trawl_b_two(param, h))
}

gamma_trawl$trawl_b_two <- function(param, h) {
  stopifnot(length(param) == 2)
  stopifnot(all(param > 0))
  a_total <- gamma_trawl$leb_a(param)
  log_scaling_factor <- (1 - param[2]) * log(1 + h / param[1])
  return(a_total * exp(log_scaling_factor))
}

gamma_trawl$trawl_b_three <- function(param, h) {
  return(gamma_trawl$trawl_b_one(param, h))
}

gamma_trawl$config <- function() {
  return(list(n_params = 2, lower = c(0.5, 2.0), upper = c(3, 4)))
}
