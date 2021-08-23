gamma_trawl <- new.env()

# param is (alpha, H)

gamma_trawl$leb_a <- function(param) {
  return(param[1] / (param[2] - 1))
}

gamma_trawl$trawl_b_one <- function(param, h) {
  return(sup_ig_trawl$leb_a(param) - gamma_trawl$trawl_b_two(param, h))
}


gamma_trawl$trawl_b_two <- function(param, h) {
  return(exp((1 - param[2]) * log(1 + h / param[1])) * gamma_trawl$leb_a(param))
}


gamma_trawl$trawl_b_three <- function(param, h) {
  return(sup_ig_trawl$trawl_b_one(param, h))
}

gamma_trawl$config <- function() {
  return(list(n_params = 2, lower = c(0.5, 2.0), upper = c(3, 4)))
}
