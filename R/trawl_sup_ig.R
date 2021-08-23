sup_ig_trawl <- new.env()

# param is (gamma, delta)

sup_ig_trawl$leb_a <- function(param) {
  return(param[1] / param[2])
}

sup_ig_trawl$trawl_b_one <- function(param, h) {
  return(
    sup_ig_trawl$leb_a(param) *
      (1 - exp(param[1] * param[2] *
        (1 - sqrt(1 + 2 * h / param[1]^2)))
      )
  )
}

sup_ig_trawl$trawl_b_two <- function(param, h) {
  return(
    exp(param[1] * param[2] *
      (1 - sqrt(1 + 2 * h / param[1]^2))) *
      sup_ig_trawl$leb_a(param)
  )
}

sup_ig_trawl$trawl_b_three <- function(param, h) {
  return(sup_ig_trawl$trawl_b_one(param, h))
}

sup_ig_trawl$config <- function() {
  return(list(n_params = 2, lower = rep(1e-5, 2), upper = rep(10.0, 2)))
}
