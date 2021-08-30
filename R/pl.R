pairwise_likelihood <- new.env()

#' @param elems Elements.
pairwise_likelihood$check_all_non_positive <- function(elems) {
  all(elems <= 0.0)
}

#' @param elems Elements.
pairwise_likelihood$check_all_positive <- function(elems) {
  all(elems > 0.0)
}

#' @param alpha Alpha parameter.
#' @param elems Elements.
pairwise_likelihood$stand_trawl_terms <- function(alpha, elems) {
  a_total <- sum(elems[1:2])
  return(-alpha * elems / a_total)
}

pairwise_likelihood$init_guess <- function(data, depth,
                                           n_trials, type = "exp") {
  cfg <- get_trawl_params_config(type)
  trawl_evaluator <- trawl_gmm$two_stage_gmm_objective(
    data = data, depth = depth, type = type
  )
  potential_param_values <- seq(
    from = cfg$lower, to = cfg$upper, length.out = n_trials
  )
  evaluator_vals <- vapply(
    potential_param_values, trawl_evaluator, 1.0
  )
  return(potential_param_values[which.min(evaluator_vals)])
}

pairwise_likelihood$pair_pdf_constructor <- function(params, type = "exp") {
  # params is (xi, sigma, kappa, trawl_params)
  b_funcs <- get_trawl_functions(type)
  b_1_func <- b_funcs[[1]]
  b_2_func <- b_funcs[[2]]
  b_3_func <- b_funcs[[3]]

  params_noven <- parametrisation_translator(
    params = params[1:3], parametrisation = "standard", target = "transform"
  )
  trawl_params <- params[4:length(params)]
  assertthat::assert_that(pairwise_likelihood$check_all_positive(trawl_params))

  return(function(xs, h) {
    jacob_cst <- 1
    return(cpp_case_separator(xs,
      alpha = params_noven[1], beta = 1.0, kappa = params_noven[3],
      b_1 = b_1_func(trawl_params, h),
      b_2 = b_2_func(trawl_params, h),
      b_3 = b_3_func(trawl_params, h)
    ) * prod(jacob_cst))
  })
}


pairwise_likelihood$parallel_apply_pl <- function(data, k, this_pl, cl) {
  n_sample <- length(data)
  xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
  return(
    sum(
      unlist(
        parallel::parApply(
          cl,
          X = xs_stack,
          MARGIN = 1,
          FUN = function(xs) {
            pl_val <- this_pl(xs, h = k)
            if (is.nan(pl_val)) {
              pl_val <- -10
            }
            if (pl_val >= 0.0) {
              pl_val <- log(pl_val + 1e-7)
            }
            return(pl_val)
          }
        )
      )
    )
  )
}

pairwise_likelihood$apply_pl <- function(data, k, this_pl) {
  n_sample <- length(data)
  xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
  return(
    sum(
      unlist(
        apply(
          X = xs_stack,
          MARGIN = 1,
          FUN = function(xs) {
            pl_val <- this_pl(xs, h = k)
            if (is.nan(pl_val)) {
              pl_val <- -10
            }
            if (pl_val >= 0.0) {
              pl_val <- log(pl_val + 1e-7)
            }
            return(pl_val)
          }
        )
      )
    )
  )
}

pairwise_likelihood$pl_constructor <- function(params, depth,
                                               pair_likehood, cl = NULL) {
  # returns function implementing Consecutive PL with depth depth
  stopifnot(depth >= 1)
  pl_f <- function(data) {
    n_sample <- length(data)
    this_pl <- compiler::cmpfun(pair_likehood)

    if (!is.null(cl)) {
      log_pl_per_depth <- vapply(seq_len(depth), # loop through depths
        FUN = function(k) {
          return(pairwise_likelihood$parallel_apply_pl(data, k, this_pl, cl))
        },
        FUN.VALUE = 1.0
      )
    } else {
      log_pl_per_depth <- vapply(seq_len(depth), # loop through depths
        FUN = function(k) {
          return(pairwise_likelihood$apply_pl(data, k, this_pl))
        },
        FUN.VALUE = 1.0
      )
    }

    # adds jacobian
    jacob <- transformation_jacobian(params_std = params[1:3])
    jacob_vals_per_depth <- vapply(
      seq_len(depth), # loop through depths
      FUN = function(k) {
        xs_stack <- c(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
        jacob_vals <- jacob(xs_stack)
        jacob_vals <- sum(log(jacob_vals[!is.na(jacob_vals)] + 1e-7))
        return(jacob_vals)
      },
      FUN.VALUE = 1.0
    )
    log_pl_per_depth <- log_pl_per_depth + jacob_vals_per_depth

    return(-sum(log_pl_per_depth)) # -1 because optim minimises by default
  }

  return(pl_f)
}


pairwise_likelihood$pl_constructor_single <- function(params, k,
                                                      pair_likehood) {
  # returns function implementing Consecutive PL with depth depth
  stopifnot(k >= 1)
  pl_f <- function(data) {
    # adds jacobian
    jacob <- transformation_jacobian(params_std = params[1:3])
    jacob_vals <- jacob(data)
    jacob_vals <- sum(log(jacob_vals[!is.na(jacob_vals)] + 1e-7))

    log_pl_per_depth <- pair_likehood(data, k) + jacob_vals
    return(-log_pl_per_depth) # -1 because optim minimises by default
  }

  return(pl_f)
}

pairwise_likelihood$trawl_pl_standard <- function(params, depth,
                                                  type = "exp", cl = NULL) {
  # param with (xi, sigma, kappa, trawl_params)
  pair_likehood_f <- pairwise_likelihood$pair_pdf_constructor(
    params = params, type = type
  ) # yields a function of (xs, h)

  return(
    pairwise_likelihood$pl_constructor(
      params = params, depth = depth, pair_likehood = pair_likehood_f, cl = cl
    )
  )
}

pairwise_likelihood$trawl_pl <- function(data, depth,
                                         type = "exp", cl = NULL) {
  return(function(params) {
    pl_functional <- pairwise_likelihood$trawl_pl_standard(
      params = params,
      depth = depth,
      type = type,
      cl = cl
    ) # returns a function of data
    return(pl_functional(data))
  })
}

pairwise_likelihood$two_stage_trawl_pl <- function(data, depth,
                                                   type = "exp", cl = NULL) {
  params_univ <- composite_marginal_mle(data)

  return(function(params) {
    pl_functional <- pairwise_likelihood$trawl_pl_standard(
      params = c(params_univ, params),
      depth = depth,
      type = type,
      cl = cl
    ) # returns a function of data
    return(pl_functional(data))
  })
}

pairwise_likelihood$trawl_pl_score <- function(params, depth,
                                               type = "exp", max_length = 100) {
  # Full Score function
  return(
    function(data) {
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]

      score_per_depth <- lapply(
        seq_len(depth),
        function(k) {
          xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
          # for each pair of data, compute score
          score_vals <- apply(
            xs_stack,
            MARGIN = 1,
            FUN = function(xs) {
              # log-PL
              log_pl <- function(par) {
                pair_pdf <- pairwise_likelihood$pair_pdf_constructor(par, type)
                pl_w_jacob <- pairwise_likelihood$pl_constructor_single(
                  par, k, pair_pdf
                )
                return(-pl_w_jacob(xs))
              }
              return(pracma::grad(log_pl, x0 = params))
            }
          )

          return(t(score_vals))
        }
      )
      return(score_per_depth)
    }
  ) # list of depth items data_length x length(params)
}

pairwise_likelihood$trawl_pl_hessian <- function(params, depth,
                                                 type = "exp",
                                                 max_length = 100) {
  # Full Score function
  return(
    function(data) {
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]

      score_per_depth <- lapply(
        seq_len(depth),
        function(k) {
          xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
          t(apply(xs_stack,
            MARGIN = 1,
            FUN = function(xs) {
              log_pl <- function(par) {
                pair_pdf <- pairwise_likelihood$pair_pdf_constructor(par, type)
                pl_w_jacob <- pairwise_likelihood$pl_constructor_single(
                  par, k, pair_pdf
                )
                return(-pl_w_jacob(xs))
              }
              return(pracma::hessian(log_pl, x0 = params))
            }
          ))
        }
      )
      return(score_per_depth)
    }
  ) # list of depth items data_length x length(params)
}

pairwise_likelihood$trawl_pl_score_partial <- function(params, depth,
                                                       type = "exp",
                                                       max_length = 100) {
  # Full Score function
  return(
    function(data) {
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]
      k <- 2

      xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
      trawl_params <- params[4:length(params)]
      model_params <- params[1:3]

      score_per_depth <- lapply(
        seq_len(depth),
        function(k) {
          apply(xs_stack,
            MARGIN = 1,
            FUN = function(xs) {
              log_pl <- function(par) {
                pair_pdf <- pairwise_likelihood$pair_pdf_constructor(
                  c(model_params, par), type
                )
                pl_w_jacob <- pairwise_likelihood$pl_constructor_single(
                  c(model_params, par), k, pair_pdf
                )
                return(-pl_w_jacob(xs))
              }
              return(pracma::grad(log_pl, x0 = trawl_params))
            }
          )
        }
      )
      return(score_per_depth)
    }
  ) # list of depth items data_length x length(params)
}

pairwise_likelihood$trawl_pl_hac <- function(data, params, depth,
                                             k = 10, type = "exp",
                                             max_length = 100) {
  lk_score <- pairwise_likelihood$trawl_pl_score(
    params, depth, type, max_length
  )
  pl_score_per_depth <- lk_score(data)

  score_acf_autocov_mat <- lapply(
    pl_score_per_depth,
    function(pl_score) {
      autocovariance_matrix(pl_score, k)
    }
  )
  pl_hac <- lapply(
    score_acf_autocov_mat, function(autocov_mat) {
      make_hac(autocov_mat)
    }
  )
  return(Reduce(`+`, pl_hac)) # sum across clusters
}

pairwise_likelihood$trawl_pl_hac_partial <- function(data, params, depth,
                                                     k = 10, type = "exp",
                                                     max_length = 100) {
  # only the trawl parameters
  lk_score <- pairwise_likelihood$trawl_pl_score_partial(
    params, depth, type, max_length
  )
  pl_score_per_depth <- lk_score(data)
  score_acf_autocov_mat <- lapply(
    pl_score_per_depth,
    function(pl_score) {
      autocovariance_matrix(pl_score, k)
    }
  )
  pl_hac <- lapply(score_acf_autocov_mat, function(autocov_mat) {
    make_hac(autocov_mat)
  })
  return(Reduce(`+`, pl_hac)) # sum across clusters
}

pairwise_likelihood$two_stage_variance <- function(data, params, depth,
                                                   type = "exp",
                                                   max_length = 100) {
  # only the trawl parameters
  lk_score <- pairwise_likelihood$trawl_pl_score(
    params, depth, type, max_length
  )
  pl_score_per_depth <- lk_score(data)

  lk_hessian <- pairwise_likelihood$trawl_pl_hessian(
    params, depth, type, max_length
  )
  pl_hessian_per_depth <- lk_hessian(data)
  mean_hessian <- Reduce(`+`, lapply(pl_hessian_per_depth, function(x) {
    apply(x, MARGIN = 2, FUN = mean)
  })) # length(data) x length(params)^2
  mean_hessian <- mean_hessian
  mean_hessian_wo_trawl <- matrix(
    mean_hessian[1:(length(params) * 3)],
    nrow = length(params), ncol = 3
  )
  mean_hessian_only_trawl <- matrix(
    mean_hessian[(length(params) * 3 + 1):length(params)^2],
    nrow = 1, ncol = length(params)
  )

  clk_hessian <- composite_likelihood_hessian(
    params[1:3],
    max_length = max_length * depth
  )
  value_clk_hessian <- clk_hessian(data)
  value_clk_hessian <- apply(value_clk_hessian, 2, mean)
  value_clk_hessian <- matrix(value_clk_hessian, 3, 3, byrow = F)

  clk_score <- composite_likelihood_score(params[1:3], max_length = max_length)
  value_clk_score <- clk_score(data)
  value_clk_score <- t(value_clk_score)

  # correct scores
  # max_length x length(params)
  correction_composite <- mean_hessian_wo_trawl %*% value_clk_hessian
  correction_composite <- correction_composite %*% value_clk_score
  correction_composite <- t(correction_composite)
  n_correction <- nrow(correction_composite)

  pair_corrections <- lapply(seq_len(depth), function(i) {
    .5 * correction_composite[1:(n_correction - i), ] +
      .5 * correction_composite[(i + 1):(n_correction), ]
  })
  corrected_per_depth <- Map("-", pl_score_per_depth, pair_corrections)
  core_with_corrections <- lapply(corrected_per_depth, function(x) {
    (t(x) %*% x) / nrow(x)
  }) # this computes E((s_pl - corr) x (s_pl - corr)^T) nolint

  core_with_corrections <- lapply(core_with_corrections, function(x) {
    mean_hessian_only_trawl %*% x %*% t(mean_hessian_only_trawl)
  })
  core_with_corrections <- mean(unlist(core_with_corrections)) # mean across k

  return(core_with_corrections)
}
