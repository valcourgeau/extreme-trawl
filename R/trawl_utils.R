get_trawl_env <- function(type) {
  return(switch(type,
    "exp" = exponential_trawl,
    "sum_exp" = sum_exponential_trawl,
    "sup_ig" = sup_ig_trawl,
    "gamma" = gamma_trawl
  ))
}

#' Obtaining the three functions necessary for the trawls.
#' @param type Trawl type (e.g. `"exp"`, `"sum_exponential_trawl"`, etc).
#' @return A vector of three functions to define trawl intersections,
#'     which depend on the trawl environment.
#' @examples
#' get_trawl_functions("exp") # Exponential trawl
get_trawl_functions <- function(type) {
  select_env <- get_trawl_env(type)

  return(
    c(
      select_env$trawl_b_one, select_env$trawl_b_two, select_env$trawl_b_three
    )
  )
}

#' Returns triplet list (n_params, lower, upper) for each trawl type.
#' @param type Trawl type (e.g. `"exp"`, `"sum_exponential_trawl"`, etc).
#' @return A list of three vectors `list(n_params, lower, upper)` to define
#'     trawl intersections parameters.
#' @examples
#' get_trawl_params_config("exp") # Exponential trawl cfg
get_trawl_params_config <- function(type) {
  return(get_trawl_env(type)$config())
}

get_trawl_envs_list <- function() {
  return(c(
    "exponential_trawl",
    "sum_exponential_trawl",
    "sup_ig_trawl",
    "gamma_trawl"
  ))
}
