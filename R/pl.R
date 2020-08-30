PairwiseLikelihood <- new.env()

PairwiseLikelihood$CheckAllNonpositive <- function(elems){all(elems <= 0.0)}

PairwiseLikelihood$CheckAllPositive <- function(elems){all(elems > 0.0)}

PairwiseLikelihood$StandTrawlTerms <- function(alpha, elems){A <- sum(elems[1:2]); return(-alpha*elems/A)}

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
      return(ev.trawl.cpp::CppCaseSeparator(xs,
                                            alpha = params_noven[1],
                                            beta = params_noven[2],
                                            kappa = params_noven[3],
                                            B1 = B1_func(trawl_params, h),
                                            B2 = B2_func(trawl_params, h),
                                            B3 = B3_func(trawl_params, h)) / abs(params[1])^3
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
              return(log(max(pl_val, 1e-9)))
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
              return(log(max(pl_val, 1e-9)))
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

PairwiseLikelihood$TrawlPLStandard <- function(params, depth, type='exp', cl=NULL){
  # param with (xi, sigma, kappa, trawl_params)
  pair_likehood_f <- PairPDFConstructor(params_noven = params_tmp, type = type) # yields a function of (xs, h)

  return(PairwiseLikelihood$PLConstructor(params = params, depth = depth, pair_likehood = pair_likehood_f, cl=cl))
}

PairwiseLikelihood$TrawlPL <- function(data, depth, type='exp', cl=NULL){
  return(function(params){
    pl_functional <- TrawlPLStandard(
      params = params,
      depth = depth,
      type=type,
      cl=cl) # returns a function of data
    return(pl_functional(data))
  })
}

PairwiseLikelihood$TwoStageTrawlPL <- function(data, depth, type='exp', cl=NULL){
  params_univ <- CustomMarginalMLE(data)

  return(function(params){
    pl_functional <- TrawlPLStandard(
      params = c(params_univ, params),
      depth = depth,
      type=type,
      cl=cl) # returns a function of data
    return(pl_functional(data))
  })
}



