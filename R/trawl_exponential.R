exponential_trawl <- new.env()

exponential_trawl$trawl_b_one <- function(param, h) {
  stopifnot(length(param) == 1)
  return((1.0 - exp(-param * h)) / param)
}

exponential_trawl$trawl_b_two <- function(param, h) {
  stopifnot(length(param) == 1)
  return(exp(-param * h) / param)
}

exponential_trawl$trawl_b_three <- function(param, h) {
  return(exponential_trawl$trawl_b_one(param, h))
}

exponential_trawl$config <- function() {
  return(list(n_params = 1, lower = 0.05, upper = 0.99))
}
