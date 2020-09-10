PairwiseLikelihood <- new.env()

PairwiseLikelihood$CheckAllNonpositive <- function(elems){all(elems <= 0.0)}

PairwiseLikelihood$CheckAllPositive <- function(elems){all(elems > 0.0)}

PairwiseLikelihood$StandTrawlTerms <- function(alpha, elems){A <- sum(elems[1:2]); return(-alpha*elems/A)}

PairwiseLikelihood$InitGuess <- function(data, depth, n_trials, type='exp'){
  cfg <- GetTrawlParamsConfig(type)
  trawl_evaluator <- TrawlGMM$TwoStageGMMObjective(data=data, depth=depth, type=type)
  potential_param_values <- seq(from=cfg$lower, to=cfg$upper, length.out=n_trials)
  evaluator_vals <- vapply(potential_param_values, trawl_evaluator, 1.0)
  return(potential_param_values[which.min(evaluator_vals)])
}

PairwiseLikelihood$PairPDFConstructor <- function(params, type='exp'){
  # params is (xi, sigma, kappa, trawl_params)
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  B3_func <- B_funcs[[3]]

  params_noven <- ParametrisationTranslator(params[1:3], parametrisation='standard', target='transform')
  trawl_params <- params[4:length(params)]
  assertthat::assert_that(PairwiseLikelihood$CheckAllPositive(trawl_params))

  return(function(xs, h){
      # are_non_zero <- abs(xs) > 1e-8
      # jacob_cst <- min(abs(params[1])^{-3}, 10000, na.rm = T) #* are_non_zero
      # jacob_cst[!are_non_zero] <- 1
      jacob_cst <- 1
      return(CppCaseSeparator(xs,
                              alpha = params_noven[1],
                              beta = 1.0,
                              kappa = params_noven[3],
                              B1 = B1_func(trawl_params, h),
                              B2 = B2_func(trawl_params, h),
                              B3 = B3_func(trawl_params, h)) * prod(jacob_cst)
      )
    }
  )
}


PairwiseLikelihood$ParallelApplyPL <- function(data, k, this_pl, cl){
  n_sample <- length(data)
  xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
  return(
    sum(
      unlist(
        parallel::parApply(
          cl,
          X = xs_stack,
          MARGIN = 1,
          FUN = function(xs){
            pl_val <- this_pl(xs, h=k)
            if(is.nan(pl_val)){cat('NA', xs, '\n'); return(-10)}
            if(pl_val < 0.0){
              return(pl_val)
            }else{
              return(log(pl_val + 1e-7))
            }})
      )
    ))
}

PairwiseLikelihood$ApplyPL <- function(data, k, this_pl){
  n_sample <- length(data)
  xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
  return(
    sum(
      unlist(
        apply(
          X = xs_stack,
          MARGIN = 1,
          FUN = function(xs){
            pl_val <- this_pl(xs, h=k)
            if(is.nan(pl_val)){cat('NA', xs, '\n'); return(-10)}
            if(pl_val < 0.0){
              return(pl_val)
            }else{
              return(log(pl_val + 1e-7))
            }})
      )
    ))
}

PairwiseLikelihood$PLConstructor <- function(params, depth, pair_likehood, cl=NULL){
  # returns function implementing Consecutive PL with depth depth
  stopifnot(depth >= 1)
  pl_f <- function(data){
    n_sample <- length(data)
    this_pl <- compiler::cmpfun(pair_likehood)

    if(!is.null(cl)){
      log_pl_per_depth <- vapply(1:depth, # loop through depths
                                 FUN = function(k){
                                    return(PairwiseLikelihood$ParallelApplyPL(data, k, this_pl, cl))
                                 },
                                 FUN.VALUE = 1.0)
    }else{
      log_pl_per_depth <- vapply(1:depth, # loop through depths
                                 FUN = function(k){
                                   return(PairwiseLikelihood$ApplyPL(data, k, this_pl))
                                 },
                                 FUN.VALUE = 1.0)
    }

    # adds jacobian
    jacob <- TransformationJacobian(params_std = params[1:3])
    jacob_vals_per_depth <-  vapply(
      1:depth, # loop through depths
      FUN = function(k){
        xs_stack <- c(data[1:(n_sample-k)], data[(k+1):(n_sample)])
        jacob_vals <- jacob(xs_stack)
        jacob_vals <- sum(log(jacob_vals[!is.na(jacob_vals)] + 1e-7))
        return(jacob_vals)
      },
      FUN.VALUE = 1.0)
    log_pl_per_depth <- log_pl_per_depth + jacob_vals_per_depth

    return(-sum(log_pl_per_depth)) # -1 because optim minimises by default
  }

  return(pl_f)
}


PairwiseLikelihood$PLConstructorSingle <- function(params, k, pair_likehood){
  # returns function implementing Consecutive PL with depth depth
  stopifnot(k >= 1)
  pl_f <- function(data){
    n_sample <- length(data)

    # adds jacobian
    jacob <- TransformationJacobian(params_std = params[1:3])
    jacob_vals <- jacob(data)
    jacob_vals <- sum(log(jacob_vals[!is.na(jacob_vals)] + 1e-7))

    log_pl_per_depth <- pair_likehood(data, k) + jacob_vals
    return(-log_pl_per_depth) # -1 because optim minimises by default
  }

  return(pl_f)
}

PairwiseLikelihood$TrawlPLStandard <- function(params, depth, type='exp', cl=NULL){
  # param with (xi, sigma, kappa, trawl_params)
  pair_likehood_f <- PairwiseLikelihood$PairPDFConstructor(params = params, type = type) # yields a function of (xs, h)

  return(PairwiseLikelihood$PLConstructor(params = params, depth = depth, pair_likehood = pair_likehood_f, cl=cl))
}

PairwiseLikelihood$TrawlPL <- function(data, depth, type='exp', cl=NULL){
  return(function(params){
    pl_functional <- PairwiseLikelihood$TrawlPLStandard(
      params = params,
      depth = depth,
      type=type,
      cl=cl) # returns a function of data
    return(pl_functional(data))
  })
}

PairwiseLikelihood$TwoStageTrawlPL <- function(data, depth, type='exp', cl=NULL){
  params_univ <- CompositeMarginalMLE(data)

  return(function(params){
    pl_functional <- PairwiseLikelihood$TrawlPLStandard(
      params = c(params_univ, params),
      depth = depth,
      type=type,
      cl=cl) # returns a function of data
    return(pl_functional(data))
  })
}

PairwiseLikelihood$TrawlPLScore <- function(params, depth, type='exp', max_length=100){
  # Full Score function
  return(
    function(data){
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]

      score_per_depth <- lapply(
         1:depth,
         function(k){
           xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
            t(apply(xs_stack, MARGIN = 1,
                   FUN = function(xs){
                     log_pl <- function(par){
                       pair_pdf <- PairwiseLikelihood$PairPDFConstructor(par, type)
                       pl_w_jacob <- PairwiseLikelihood$PLConstructorSingle(par, k, pair_pdf)
                       return(-pl_w_jacob(xs))
                     }
                     return(pracma::grad(log_pl, x0 = params))}
                   )
              )
         }
      )
      return(score_per_depth)
    }) # list of depth items data_length x length(params)
}

PairwiseLikelihood$TrawlPLHessian <- function(params, depth, type='exp', max_length=100){
  # Full Score function
  return(
    function(data){
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]

      score_per_depth <- lapply(
        1:depth,
        function(k){
          xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
          t(apply(xs_stack, MARGIN = 1,
                  FUN = function(xs){
                    log_pl <- function(par){
                      pair_pdf <- PairwiseLikelihood$PairPDFConstructor(par, type)
                      pl_w_jacob <- PairwiseLikelihood$PLConstructorSingle(par, k, pair_pdf)
                      return(-pl_w_jacob(xs))
                    }
                    return(pracma::hessian(log_pl, x0 = params))}
          )
          )
        }
      )
      return(score_per_depth)
    }) # list of depth items data_length x length(params)
}

PairwiseLikelihood$TrawlPLScorePartial <- function(params, depth, type='exp', max_length=100){
  # Full Score function
  return(
    function(data){
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]
      k <- 2

      xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
      trawl_params <- params[4:length(params)]
      model_params <- params[1:3]

      score_per_depth <- lapply(
        1:depth,
        function(k){
          apply(xs_stack, MARGIN = 1,
                  FUN = function(xs){
                    log_pl <- function(par){
                      pair_pdf <- PairwiseLikelihood$PairPDFConstructor(c(model_params, par), type)
                      pl_w_jacob <- PairwiseLikelihood$PLConstructorSingle(c(model_params, par), k, pair_pdf)
                      return(-pl_w_jacob(xs))
                    }
                    return(pracma::grad(log_pl, x0 = trawl_params))}
          )
        }
      )
      return(score_per_depth)
    }) # list of depth items data_length x length(params)
}

PairwiseLikelihood$TrawlPLHAC <- function(data, params, depth, k=10, type='exp', max_length=100){
  lk_score <- PairwiseLikelihood$TrawlPLScore(params, depth, type, max_length)
  pl_score_per_depth <- lk_score(data)

  score_acf_autocov_mat <- lapply(
    pl_score_per_depth,
    function(pl_score){AutocovarianceMatrix(pl_score, params, k)})
  pl_hac <- lapply(score_acf_autocov_mat, function(autocov_mat){MakeHAC(autocov_mat)})
  return(Reduce(`+`, pl_hac)) # sum across clusters
}

PairwiseLikelihood$TrawlPLHACPartial <- function(data, params, depth, k=10, type='exp', max_length=100){
  # only the trawl parameters
  lk_score <- PairwiseLikelihood$TrawlPLScorePartial(params, depth, type, max_length)
  pl_score_per_depth <- lk_score(data)

  trawl_params <- params[4:length(params)]

  score_acf_autocov_mat <- lapply(
    pl_score_per_depth,
    function(pl_score){AutocovarianceMatrix(pl_score, trawl_params, k)})
  pl_hac <- lapply(score_acf_autocov_mat, function(autocov_mat){MakeHAC(autocov_mat)})
  return(Reduce(`+`, pl_hac)) # sum across clusters
}

PairwiseLikelihood$TwoStageVariance <- function(data, params, depth, type='exp', max_length=100){
  # only the trawl parameters
  lk_score <- PairwiseLikelihood$TrawlPLScore(params, depth, type, max_length)
  pl_score_per_depth <- lk_score(data)

  lk_hessian <- PairwiseLikelihood$TrawlPLHessian(params, depth, type, max_length)
  pl_hessian_per_depth <- lk_hessian(data)
  mean_hessian <- Reduce(`+`, lapply(pl_hessian_per_depth, function(x){apply(x, MARGIN = 2, FUN = mean)})) # length(data) x length(params)^2
  mean_hessian <- mean_hessian
  mean_hessian_wo_trawl <- matrix(mean_hessian[1:(length(params)*3)], nrow = length(params), ncol=3)
  mean_hessian_only_trawl <- matrix(mean_hessian[(length(params)*3+1):length(params)^2], nrow = 1, ncol=length(params))

  print(mean_hessian_wo_trawl)
  print(mean_hessian_only_trawl)

  clk_hessian <- CompositeLikelihoodHessian(params[1:3], max_length=max_length*depth)
  value_clk_hessian <- clk_hessian(data)
  value_clk_hessian <- apply(value_clk_hessian, 2, mean)
  value_clk_hessian <- matrix(value_clk_hessian, 3, 3, byrow = F)
  # value_clk_hessian <- solve(value_clk_hessian)

  clk_score <- CompositeLikelihoodScore(params[1:3], max_length=max_length)
  value_clk_score <- clk_score(data)
  value_clk_score <- t(value_clk_score)

  # correct scores
  correction_composite <- t(mean_hessian_wo_trawl %*% value_clk_hessian %*% value_clk_score) # max_length x length(params)
  n_correction <- nrow(correction_composite)
  pair_corrections <- lapply(1:depth, function(i){.5*correction_composite[1:(n_correction-i),] + .5*correction_composite[(i+1):(n_correction),]})
  corrected_per_depth <- Map('-', pl_score_per_depth, pair_corrections)
  core_with_corrections <- lapply(corrected_per_depth, function(x){(t(x) %*% x) / nrow(x)}) # E((s_pl - corr) * (s_pl - corr)^T)

  core_with_corrections <- lapply(core_with_corrections, function(x){mean_hessian_only_trawl %*% x %*% t(mean_hessian_only_trawl)})
  core_with_corrections <- mean(unlist(core_with_corrections)) # mean across k

  return(core_with_corrections)
}
