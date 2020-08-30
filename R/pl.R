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

  params_noven <- ParametrisationTranslator(params, parametrisation='standard', target='transform')

  trawl_params <- params_noven[4:length(params_noven)]
  assertthat::assert_that(CheckAllPositive(trawl_params))
  return(function(xs, h){
    return(ev.trawl.cpp::CppCaseSeparator(xs,
                                          alpha = params_noven[1],
                                          beta = params_noven[2],
                                          kappa = params_noven[3],
                                          B1 = B1_func(trawl_params, h),
                                          B2 = B2_func(trawl_params, h),
                                          B3 = B3_func(trawl_params, h) / abs(params[1])^3
    ))
  })
}

PairwiseLikelihood$ParallelApplyPL <- function(data, k, this_pl){
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


PairwiseLikelihood$PLConstructor <- function(params, depth, pair_likehood, cl=NULL){
  # returns function implementing Consecutive PL with depth depth
  stopifnot(depth >= 1)
  pl_f <- function(data){
    n_sample <- length(data)
    this_pl <- cmpfun(pair_likehood)

    if(!is.null(cl)){
      log_pl_per_depth <- vapply(1:depth, # loop through depths
                                 FUN = function(k){
                                    return(PairwiseLikelihood$ParallelApplyPL(data, k, this_pl))
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
    jacob <- TransformationJacobian(params_std = params)
    jacob_vals_per_depth <-  vapply(
      1:depth, # loop through depths
      FUN = function(k){
        xs_stack <- c(data[1:(n_sample-k)], data[(k+1):(n_sample)])
        jacob_vals <- jacob(data_for_mle)
        jacob_vals <- sum(log(jacob_vals[!is.na(jacob_vals)] + 1e-7))
        return(jacob_vals)
      },
      FUN.VALUE = 1.0)
    log_pl_per_depth <- log_pl_per_depth + jacob_vals_per_depth

    return(sum(log_pl_per_depth))
  }

  return(pl_f)
}

PairwiseLikelihood$TrawlPLStandard <- function(params, depth, type='exp', cl=NULL){
  # param with (xi, sigma, kappa, trawl_params)
  B_funcs <- GetTrawlFunctions(type)
  B1_func <- B_funcs[[1]]
  B2_func <- B_funcs[[2]]
  B3_func <- B_funcs[[3]]

  pair_likehood_f <- PairPDFConstructor(params_noven = params_tmp, type = type) # yields a function of (xs, h)

  return(PairwiseLikelihood$PLConstructor(params = params, depth = depth, pair_likehood = pair_likehood_f, cl=cl))
}

PairwiseLikelihood$TrawlPLFunctional <- function(params, depth, type='exp', cl=NULL){
  PLOperator <- PairwiseLikelihood$TrawlPLStandard(params = params, depth = depth, type = type, cl=cl)

  # multiply by - 1
  return(function(data){
    return((-1)*PLOperator(data))})
}

PairwiseLikelihood$TrawlPL <- function(data, depth, type='exp', parametrisation='standard', cl=NULL){
  return(function(params){
    pl_functional <- TrawlPLFunctional(params = params,
                                       depth = depth,
                                       type=type,
                                       parametrisation=parametrisation,
                                       cl=cl) # returns a function of data
    return(pl_functional(data))
  })
}



