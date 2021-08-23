
grid_foundations <- function(n, vanishing_depth, x = 1) {
  # returns a triangular matrix
  n <- n + vanishing_depth # accommodate for warm-up

  n_elems_per_line <- vanishing_depth
  index_i <- unlist(lapply(1:max(2, n_elems_per_line), FUN = function(i) {
    return(1:(n))
  }))
  index_j <- unlist(lapply(1:max(2, n_elems_per_line), FUN = function(i) {
    return(i:(n + i - 1))
  }))
  return(Matrix::sparseMatrix(i = index_i, j = index_j, use.last.ij = T, x = x))
}

trawl_slicing <- function(n, vanishing_depth, trawl_parameter, type = "exp") {
  # returns matrix with
  # [[S(1,1), S(2,1), S(3,1), ..., S(vanishing_depth,1)],
  #  [S(2,2), S(3,2), S(4,2), ..., S(vanishing_depth+1,2)],
  #  ...,
  #  [S(n,n), S(3,2), S(4,2), ..., S(vanishing_depth+n,n)],
  # ]

  b_funcs <- get_trawl_functions(type = type)
  b_2_func <- b_funcs[[2]]

  one_split <- diff(-b_2_func(param = trawl_parameter, h = 0:vanishing_depth))
  slices <- matrix(rep(one_split, n), ncol = vanishing_depth, byrow = T)
  a_total <- b_2_func(trawl_parameter, 0.0)
  slices[2:nrow(slices), ] <- slices[2:nrow(slices), ] - rep(1, n - 1) %o%
    c(one_split[2:length(one_split)], 0.0)

  return(slices / a_total) # divide by \mu^{leb}(A)
}


gamma_grid <- function(alpha, beta, n,
                       vanishing_depth, trawl_parameter, type = "exp") {
  n_diags <- max(2, vanishing_depth)
  n_non_zero_elements <- n_diags * n

  gamma_shapes <- generate_shapes(
    alpha = alpha,
    n = n,
    vanishing_depth = vanishing_depth,
    trawl_parameter = trawl_parameter,
    type = type
  )
  stopifnot(length(gamma_shapes) == n_non_zero_elements)
  gamma_sim_vals <- rgamma(
    n = n_non_zero_elements,
    shape = gamma_shapes,
    rate = beta
  )
  return(grid_foundations(
    n = n, vanishing_depth = vanishing_depth, x = gamma_sim_vals
  ))
}

generate_shapes <- function(alpha, n, vanishing_depth,
                            trawl_parameter, type = "exp") {
  stopifnot(length(alpha) == 1)

  slices <- trawl_slicing(
    n = n, vanishing_depth = vanishing_depth,
    trawl_parameter = trawl_parameter, type = type
  ) # already standardised
  return(alpha * as.vector(slices))
}

block_index <- function(n_block, n, vanishing_depth) {
  return(
    list(
      block_i = max(c(1, n_block - vanishing_depth + 1)):min(c(n, n_block)),
      block_j = n_block:min(c(n, vanishing_depth + n_block - 1))
    )
  )
}

gamma_orchestra <- function(scaled_gamma_grid, parallel = T) {
  n <- dim(scaled_gamma_grid)[1]
  vanishing_depth <- dim(scaled_gamma_grid)[2] - n + 1

  if (parallel) {
    cores <- parallel::detectCores(logical = TRUE)
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("block_index", "n", "vanishing_depth"))
    parallel::clusterEvalQ(cl, library(Matrix))
    tmp <- parallel::parLapply(
      cl = cl,
      X = n:1,
      fun = function(i) {
        blck_ind <- block_index(i, n = n, vanishing_depth = vanishing_depth)
        return(sum(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j]))
      }
    )
    parallel::stopCluster(cl)
    return(unlist(tmp))
  } else {
    return(vapply(1:n, FUN = function(i) {
      blck_ind <- block_index(i, n = n, vanishing_depth = vanishing_depth)
      return(sum(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j]))
    }, FUN.VALUE = 1.0))
  }
}

#' Returns the coverage of the trawl set approximation.
#' @param trawl_parameter Trawl parameter.
#' @param vanishing_depth Approximation depth.
#' @param type Trawl type (e.g. `"exp"`, `"sum_exp"`).
#' @param get_value Prints or returns the value. Defaults to `FALSE`.
#' @return Trawl set approximation (print or value).
#' @examples
#' n <- 10
#' vanishing_depth <- 3
#' gf <- grid_foundations(n, vanishing_depth)
#'
#' tmp <- gamma_grid(
#'   alpha = 3, beta = 20, n = 2500,
#'   vanishing_depth = 30, trawl_parameter = 0.1
#' )
#' tmp
#' plot(gamma_orchestra(tmp))
#' b_funcs <- get_trawl_functions(type = type)
#' b_1_func <- b_funcs[[1]]
#' b_2_func <- b_funcs[[2]]
#' b_3_func <- b_funcs[[3]]
#' b_2_func(0.3, 0:10)
#' (-b_2_func(param = 0.8, h = 0:k)) %>% diff()
#' print_vanishing_coverage(trawl_parameter = 0.3, vanishing_depth = 10)
print_vanishing_coverage <- function(trawl_parameter, vanishing_depth,
                                     type = "exp", get_value = F) {
  b_funcs <- get_trawl_functions(type = type)
  b_2_func <- b_funcs[[2]]

  max_val <- b_2_func(trawl_parameter, h = 0.0)

  coverage <- -sum(
    diff(
      abs(
        (b_2_func(param = trawl_parameter, h = c(0, vanishing_depth)))
      )
    )
  ) / abs(max_val)
  if (get_value) {
    return(coverage)
  } else {
    cat("Coverage:", round(coverage * 100.0, 2), "%\n")
  }
}


trawl_simulation <- function(alpha, beta, n, vanishing_depth,
                             trawl_parameter, type, parallel = F) {
  gamma_grid <- gamma_grid(
    alpha = alpha,
    beta = beta,
    n = n,
    vanishing_depth = vanishing_depth,
    trawl_parameter = trawl_parameter,
    type = type
  )
  return(gamma_orchestra(gamma_grid, parallel = parallel))
}

# nolint start
exceedances_simulation <- function(params, n, vanishing_depth, type,
                                   m = 1, parametrisation = "standard",
                                   parallel = F, algo = "standard") {
  params_trawl <- parametrisation_translator(
    params = params[1:3], parametrisation = parametrisation,
    target = "transform"
  )
  kappa <- params[3]
  trawl_parameter <- params[4:length(params)]

  coverage <- print_vanishing_coverage(
    trawl_parameter = trawl_parameter, vanishing_depth = vanishing_depth,
    type = type, get_value = T
  )

  if (parallel) {
    cores <- parallel::detectCores(logical = TRUE)
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(
      cl, c(
        "transformation_map_inverse",
        "transformation_map",
        "transformation_jacobian",
        "parametrisation_translator",
        "pairwise_likelihood",
        "composite_marginal_mle",
        "trawl_gmm",
        "trawl_autocorrelation",
        "ev_trawl_fit",
        "trawl_simulation",
        "print_vanishing_coverage",
        "gamma_orchestra",
        "block_index",
        "grid_foundations",
        "trawl_slicing",
        "gamma_grid",
        "generate_shapes",
        "exceedances_simulation",
        get_trawl_envs_list()
      )
    )

    sim_fn <- function(i) {
      exceedances_simulation(
        params = params, n = n,
        vanishing_depth = vanishing_depth, type = type,
        m = 1, parametrisation = parametrisation,
        parallel = F, algo = algo
      )
    }

    parallel::clusterExport(
      cl, c("n", "m", "vanishing_depth", "type")
    )


    sims <- parallel::parLapply(
      cl = cl,
      X = 1:m,
      fun = sim_fn
    )
    parallel::stopCluster(cl)
    return(sims)
  }

  if (algo == "standard") {
    print("standard algo")
    trawl_simulation <- trawl_simulation(
      alpha = 1.,
      beta = 1.,
      trawl_parameter = trawl_parameter,
      n = n,
      vanishing_depth = vanishing_depth,
      type = type,
      parallel = parallel
    )
    probabilities_zero <- 1 - exp(-kappa * trawl_simulation)
    uniform_samples <- runif(n = n, min = 0, max = 1)
    who_is_extreme <- uniform_samples > probabilities_zero

    exceedances <- rep(0, n)
    exceedances[who_is_extreme] <- rexp(
      n = length(which(who_is_extreme)),
      rate = trawl_simulation[who_is_extreme]
    )
    exceedances[who_is_extreme] <- transformation_map(
      x = exceedances[who_is_extreme],
      params_std = params[1:3]
    )

    return(list(
      exceedances = exceedances,
      latent = trawl_simulation,
      coverage = coverage
    ))
  } else {
    if (algo == "cross") {
      print("cross algo")


      sim_data <- lapply(
        1:m,
        function(i) {
          exceedances_simulation(
            params = params, parametrisation = parametrisation,
            n = n, m = m, vanishing_depth = vanishing_depth,
            type = type, parallel = F, algo = "standard"
          )
        }
      )

      backbone_data <- exceedances_simulation(
        params = params, parametrisation = parametrisation,
        n = n, m = m, vanishing_depth = vanishing_depth,
        type = type, parallel = F, algo = "standard"
      )
      backbone_exc <- backbone_data$exceedances

      sim_data <- abind::abind(lapply(sim_data, function(x) {
        do.call(cbind, x)
      }), along = 3) # [,2,] for latent; [,1,] for exceedances
      sim_exc <- sim_data[, 1, ]

      exceedances <- rep(0, n)
      print(apply(sim_exc[which(backbone_exc > 0), ], 1, function(x) {
        mean(x[x > 0])
      }))
      exceedances[which(backbone_exc > 0)] <- unlist(
        apply(
          sim_exc[which(backbone_exc > 0), ], 1, function(x) {
            mean(x[x > 0])
          }
        )
      )

      acf(exceedances)
      return(list(exceedances = exceedances, coverage = coverage))
    }
    if (algo == "corr_unif") {
      trawl_simulation <- trawl_simulation(
        alpha = 1.,
        beta = 1.,
        trawl_parameter = trawl_parameter,
        n = n,
        vanishing_depth = vanishing_depth,
        type = type,
        parallel = parallel
      )
      trawl_simulation_unif <- trawl_simulation(
        alpha = 1.,
        beta = 1.,
        trawl_parameter = trawl_parameter,
        n = n,
        vanishing_depth = vanishing_depth,
        type = type,
        parallel = parallel
      )

      probabilities_zero <- 1 - exp(-kappa * trawl_simulation)
      corr_uniform <- pgamma(trawl_simulation_unif, shape = 1, rate = 1)

      who_is_extreme <- corr_uniform > probabilities_zero

      exceedances <- rep(0, n)
      exceedances[who_is_extreme] <- rexp(
        n = length(which(who_is_extreme)),
        rate = trawl_simulation[who_is_extreme]
      )
      exceedances[who_is_extreme] <- transformation_map(
        x = exceedances[who_is_extreme],
        params_std = params[1:3]
      )

      return(list(
        exceedances = exceedances,
        latent = trawl_simulation,
        coverage = coverage
      ))
    }
    if (algo == "dynamic_latent" | algo == "dynamic_uniform") {
      trawl_simulation <- trawl_simulation(
        alpha = 1.,
        beta = 1.,
        trawl_parameter = trawl_parameter,
        n = n,
        vanishing_depth = vanishing_depth,
        type = type,
        parallel = parallel
      )

      b_funcs <- get_trawl_functions(type = type)
      b_1_func <- b_funcs[[1]]
      b_2_func <- b_funcs[[2]]
      b_3_func <- b_funcs[[3]]

      a_value <- b_1_func(param = trawl_parameter, h = 1)
      a_value <- a_value + b_2_func(param = trawl_parameter, h = 1)
      b_1_minus_0 <- -b_1_func(param = trawl_parameter, h = 1) / a_value
      b_0_minus_1 <- -b_3_func(param = trawl_parameter, h = 1) / a_value
      b_0_1 <- -b_2_func(param = trawl_parameter, h = 1) / a_value

      probabilities_zero <- 1 - exp(-kappa * trawl_simulation)
      uniform_samples <- runif(n = n, min = 0, max = 1)

      prev_sample <- NULL
      exceedances <- apply(
        cbind(probabilities_zero, uniform_samples, trawl_simulation),
        MARGIN = 1,
        FUN = function(p_u_t) {
          if (is.null(prev_sample)) {
            extreme_proba <- 1 / (1 + kappa)
          } else {
            tmp1 <- (1 + kappa + prev_sample)^{
              1 + b_0_minus_1
            }
            tmp2 <- (1 + 2 * kappa + prev_sample)^b_0_1
            tmp3 <- (1 + kappa)^b_1_minus_0
            extreme_proba <- tmp1 * tmp2 * tmp3

            if (algo == "dynamic_latent") {
              if (prev_sample == 0.0) {
                extreme_proba <- 1 / (1 + kappa) - extreme_proba / (1 + kappa)
                extreme_proba <- extreme_proba / (1 - 1 / (1 + kappa))
              }
            } else {
              if (algo == "dynamic_uniform") {
                if (prev_sample == 0.0) {
                  extreme_proba <- extreme_proba / (1 + kappa)
                  extreme_proba <- 1 / (1 + kappa) - extreme_proba
                  extreme_proba <- extreme_proba / (1 - 1 / (1 + kappa))
                }
              }
            }
          }
          if (algo == "dynamic_latent") {
            if (1 - p_u_t[1] > extreme_proba) {
              prev_sample <<- rexp(n = 1, rate = p_u_t[3])
            } else {
              prev_sample <<- 0.0
            }
          } else {
            if (algo == "dynamic_uniform") {
              if (p_u_t[2] > extreme_proba) {
                prev_sample <<- 0.0
              } else {
                prev_sample <<- rexp(n = 1, rate = p_u_t[3])
              }
            }
          }


          return(prev_sample)
        }
      )

      who_is_extreme <- exceedances > 0

      exceedances[who_is_extreme] <- transformation_map(
        x = exceedances[who_is_extreme],
        params_std = params[1:3]
      )

      return(list(
        exceedances = exceedances,
        latent = trawl_simulation,
        coverage = coverage
      ))
    }
  }

  stop("Wrong simulation algo.")
}
# nolint end
