trawl_gmm <- new.env()

trawl_gmm$diff_fn <- function(x, y) {
  return(sum(x - y))
}

trawl_gmm$trawl_objective <- function(data, depth,
                                      type = "exp", metric = measures::SSE) {
  function(pars) {
    # pars should be (xi, sigma, kappa, rho)
    noven_pars <- parametrisation_translator(
      params = pars[1:3], parametrisation = "standard", target = "transform"
    )
    sample_cross_mom <- acf(data, plot = F, lag.max = depth)$acf

    return(function(trawl_params) {
      acf_vals <- trawl_autocorrelation$acf_trawl_collection(
        h = c(0.01, 1:(depth)), alpha = 1.,
        beta = 1., kappa = noven_pars[3],
        rho = trawl_params, delta = 0.1, cov = F, type = type
      )

      return(metric(acf_vals, sample_cross_mom))
    })
  }
}

trawl_gmm$full_gmm_objective <- function(data, depth,
                                         omega = "id", type = "exp") {
  composite <- composite_likelihood(data = data)
  trawl_objective <- trawl_gmm$trawl_objective(
    data = data,
    depth = depth,
    type = type,
    metric = trawl_gmm$diff_fn
  ) # without derivative on the last one
  return(function(par) {
    grad_vec <- c(
      pracma::grad(composite, x0 = par[1:3]) / length(data),
      trawl_objective(par)(par[4:length(par)])
    )

    if (omega == "id") {
      omega <- diag(rep(1, length(par)))
    } else {
      if (omega == "centered") {
        omega <- diag(rep(1, length(par)))
      }
    }

    return(t(grad_vec) %*% omega %*% grad_vec)
  })
}

trawl_gmm$two_stage_gmm_objective <- function(data, depth,
                                              type = "exp",
                                              metric = measures::SSE) {
  pars <- composite_marginal_mle(data = data)
  # return a function of the whole set of params
  trawl_obj <- trawl_gmm$trawl_objective(
    data, depth,
    type = type, metric = metric
  )
  return(trawl_obj(pars)) # function of trawl parameters
}

trawl_gmm$objective_single <- function(xs, k, mean_data,
                                       type = "exp",
                                       metric = measures::SSE) {
  # metric can be measures::SSE or measures::MAE
  # return the trawl objective for each k and data points

  return(
    function(pars) {
      # pars should be (xi, sigma, kappa, rho)
      noven_pars <- parametrisation_translator(
        params = pars[1:3], parametrisation = "standard", target = "transform"
      )

      return(function(trawl_params) {
        acf_vals <- trawl_autocorrelation$acf_trawl_collection(
          h = k, alpha = 1.,
          beta = 1., kappa = noven_pars[3],
          rho = trawl_params, delta = 0.1, cov = F, type = type
        )

        xs_stack <- prod(xs) - mean_data^2

        return(metric(acf_vals, xs_stack))
      })
    }
  )
}

trawl_gmm$composite_objective_single <- function(xs, k,
                                                 mean_data,
                                                 length_data,
                                                 omega = "id",
                                                 type = "exp") {
  composite <- composite_likelihood(data = xs)
  trawl_objective <- trawl_gmm$objective_single(
    xs = xs,
    k = k,
    mean_data = mean_data,
    type = type,
    metric = trawl_gmm$diff_fn
  )
  return(function(par) {
    grad_vec <- c(
      -pracma::grad(composite, x0 = par[1:3]),
      trawl_objective(par)(par[4:length(par)])
    )

    return(grad_vec)
  })
}

trawl_gmm$objective_single_hessian <- function(xs, k, mean_data, length_data,
                                               omega = "id", type = "exp") {
  composite <- composite_likelihood(data = xs)
  trawl_objective <- trawl_gmm$objective_single(
    xs = xs, k = k, mean_data = mean_data, type = type,
    metric = trawl_gmm$diff_fn
  )
  return(function(par) {
    trawl_obj <- function(p) {
      trawl_objective(c(par[1:2], p))(p[2:length(p)])
    } # p should only be (kappa, rho)
    trawl_grad <- pracma::grad(trawl_obj, x0 = par[3:length(par)])
    res <- matrix(0, length(par), length(par))
    res[4:length(par), 3:length(par)] <- trawl_grad
    res[3, 4:length(par)] <- trawl_grad[1] # grad only wrt (kappa, rho)

    composite_hessian <- -pracma::hessian(composite, x0 = par[1:3])
    res[1:3, 1:3] <- res[1:3, 1:3] + composite_hessian
    return(res)
  })
}

trawl_gmm$wrapper_objective_single <- function(xs, k, mean_data, length_data,
                                               omega = "id", type = "exp") {
  trawl_objective <- trawl_gmm$objective_single(
    xs = xs,
    k = k,
    mean_data = mean_data,
    type = type,
    metric = trawl_gmm$diff_fn
  )
  return(function(par) {
    grad_vec <- trawl_objective(par)(par[4:length(par)])
    return(grad_vec)
  })
}

trawl_gmm$trawl_gmm_score <- function(params, depth,
                                      type = "exp", max_length = 100) {
  # Full Score function
  return(
    function(data) {
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]
      mean_data <- mean(data)

      score_per_depth <- lapply(
        seq_len(depth),
        function(k) {
          xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
          t(apply(xs_stack,
            MARGIN = 1,
            FUN = function(xs) {
              full_gmm <- trawl_gmm$composite_objective_single(
                xs, k, mean_data, n_sample
              )
              return(full_gmm(params))
            }
          ))
        }
      )
      return(score_per_depth)
    }
  ) # list of depth items data_length x length(params)
}

trawl_gmm$trawl_gmm_hessian <- function(params, depth,
                                        type = "exp", max_length = 100) {
  # Full Score function
  return(
    function(data) {
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]
      mean_data <- mean(data)

      score_per_depth <- lapply(
        seq_len(depth),
        function(k) {
          xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
          t(apply(xs_stack,
            MARGIN = 1,
            FUN = function(xs) {
              hessian_gmm <- trawl_gmm$objective_single_hessian(
                xs, k, mean_data, n_sample
              )
              return(hessian_gmm(params))
            }
          ))
        }
      )
      return(score_per_depth)
    }
  ) # list of depth items data_length x length(params)
}

trawl_gmm$trawl_gmm_score_partial <- function(params, depth,
                                              type = "exp", max_length = 100) {
  # Full Score function
  return(
    function(data) {
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]
      mean_data <- mean(data)

      score_per_depth <- lapply(
        seq_len(depth),
        function(k) {
          xs_stack <- cbind(data[1:(n_sample - k)], data[(k + 1):(n_sample)])
          apply(xs_stack,
            MARGIN = 1,
            FUN = function(xs) {
              full_gmm <- trawl_gmm$wrapper_objective_single(
                xs, k, mean_data, n_sample
              )
              return(full_gmm(params))
            }
          )
        }
      )
      return(score_per_depth)
    }
  ) # list of depth items data_length x length(params)
}

trawl_gmm$trawl_gmm_hac <- function(data, params, depth,
                                    k = 10, type = "exp", max_length = 100) {
  lk_score <- trawl_gmm$trawl_gmm_score(params, depth, type, max_length)
  pl_score_per_depth <- lk_score(data)

  score_acf_autocov_mat <- lapply(
    pl_score_per_depth, function(pl_score) {
      autocovariance_matrix(pl_score, k)
    }
  )
  pl_hac <- lapply(score_acf_autocov_mat, function(autocov_mat) {
    make_hac(autocov_mat, near_pd = F)
  })
  # sum across clusters
  return(as.matrix(Matrix::nearPD(Reduce(`+`, pl_hac))$mat)) # nolint
}

trawl_gmm$trawl_gmm_hac_partial <- function(data, params, depth, k = 10,
                                            type = "exp", max_length = 100) {
  # only the trawl parameters
  lk_score <- trawl_gmm$trawl_gmm_score_partial(params, depth, type, max_length)
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


trawl_gmm$two_stage_variance <- function(data, params, depth,
                                         type = "exp", max_length = 100) {
  # return var not scaled by number of data points

  # only the trawl parameters
  lk_score <- trawl_gmm$trawl_gmm_score(params, depth, type, max_length)
  pl_score_per_depth <- lk_score(data)
  pl_score_per_depth <- lapply(
    pl_score_per_depth,
    function(x) {
      y <- x
      y[, 4:(length(params))] <- y[, 4:length(params)] / depth
      return(y)
    }
  ) # TODO RM?

  lk_hessian <- trawl_gmm$trawl_gmm_hessian(params, depth, type, max_length)
  pl_hessian_per_depth <- lk_hessian(data)
  pl_hessian_per_depth <- lapply(
    pl_hessian_per_depth,
    function(x) {
      y <- x
      y[, 4:length(params)] <- y[, 4:length(params)] / depth
      return(y)
    }
  )
  mean_hessian <- Reduce(`+`, lapply(pl_hessian_per_depth, function(x) {
    apply(x, MARGIN = 2, FUN = mean)
  }))
  mean_hessian <- mean_hessian
  mean_hessian_wo_trawl <- matrix(
    mean_hessian[1:(length(params) * 3)],
    nrow = length(params), ncol = 3
  )
  mean_hessian_only_trawl <- matrix(
    mean_hessian[(length(params) * 3 + 1):length(params)^2],
    nrow = 1, ncol = length(params)
  )

  mean_hessian_only_trawl[4:length(params)] <- mean_hessian_only_trawl[
    4:length(params)
  ] / depth

  clk_hessian <- composite_likelihood_hessian(
    params = params[1:3], max_length = max_length * depth
  )
  value_clk_hessian <- clk_hessian(data)
  value_clk_hessian <- apply(value_clk_hessian, 2, mean)
  value_clk_hessian <- matrix(value_clk_hessian, 3, 3, byrow = F)

  clk_score <- composite_likelihood_score(params[1:3], max_length = max_length)
  value_clk_score <- clk_score(data)
  value_clk_score <- t(value_clk_score)

  # correct scores: max_length x length(params)
  correction_composite <- mean_hessian_wo_trawl %*% value_clk_hessian
  correction_composite <- t(correction_composite %*% value_clk_score)

  n_correction <- nrow(correction_composite)
  pair_corrections <- lapply(seq_len(depth), function(i) {
    .5 * correction_composite[1:(n_correction - i), ] +
      .5 * correction_composite[(i + 1):(n_correction), ]
  })
  corrected_per_depth <- Map("-", pl_score_per_depth, pair_corrections)
  core_with_corrections <- lapply(corrected_per_depth, function(x) {
    (t(x) %*% x) / nrow(x)
  }) # creates E((s_pl - corr) * (s_pl - corr)^T)

  core_with_corrections <- lapply(core_with_corrections, function(x) {
    mean_hessian_only_trawl %*% x %*% t(mean_hessian_only_trawl)
  })
  core_with_corrections <- mean(unlist(core_with_corrections)) # mean across k

  return(core_with_corrections)
}
