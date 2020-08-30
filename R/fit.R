
EVTrawlFit <- function(data, depth, method, mode='two-stage', type='exp', bounds='config', cl=NULL, ...){
  # method 'PL' or 'GMM' either 'full' or 'two-stage' modes

  trawl_cfg <- GetTrawlParamsConfig(type)
  init_guess_lower_upper <- GetInitialGuessAndBounds(data)
  init_guess_model <- init_guess_lower_upper$init_guess
  lower_model <- init_guess_lower_upper$lower
  upper_model <- init_guess_lower_upper$upper

  init_trawl <- vapply(1:trawl_cfg$n_params, function(i){
    runif(n = 1,
          min = trawl_cfg$lower[i],
          max = trawl_cfg$upper[i])
  }, 1.)

  init_trawl <- vapply(1:trawl_cfg$n_params, function(i){
    runif(n = 1,
          min = trawl_cfg$lower[i],
          max = trawl_cfg$upper[i])
  }, 1.)
  if(mode == 'two-stage'){
    marginal_params <- CompositeMarginalMLE(data)
    # choose function
    if(method == 'PL'){
      optim_fn <- PairwiseLikelihood$TwoStageTrawlPL(
        data = data, depth = depth,
        type = type, cl = cl)
    }else{
      if(method == 'GMM'){
        optim_fn <- TwoStageGMMObjective(data = data, depth = depth)
      }
    }

    # choose bounds
    lower <- trawl_cfg$lower
    upper <- trawl_cfg$upper
    if(bounds == 'multipler'){
      lower <- lower * 0.8
      upper <- upper * 1.2
    }
    init_guess <- init_trawl
  }else{
    if(mode == 'full'){
      # choose function
      if(method == 'PL'){
        optim_fn <- PairwiseLikelihood$TrawlPL(
          data = data, depth = depth,
          type = type, cl = cl)
      }else{
        if(method == 'GMM'){
          optim_fn <- FullGMMObjective(
            data = data, depth = depth)
        }
      }

      # choose bounds
      lower <- trawl_cfg$lower
      upper <- trawl_cfg$upper
      if(bounds == 'multipler'){
        init_trawl <- PairwiseLikelihood$InitGuess(data=data, depth=depth, n_trials=20, type=type)
        lower <- init_trawl * 0.8
        upper <- init_trawl * 1.2
      }
      lower <- c(lower_model, lower)
      upper <- c(upper_model, upper)

      # start values
      init_guess <- c(init_guess_model, init_trawl)
    }else{
      stop(paste('mode should be two-stage or full, received', mode))
    }
  }


  trawl_inference <- stats::optim(
    fn = optim_fn,
    par = init_guess,
    lower = lower,
    upper = upper,
    method = 'L-BFGS-B',
    ...)
  if(mode == 'two-stage'){
    return(c(marginal_params, trawl_inference$par))
  }else{
    if(mode == 'full'){
      return(trawl_inference$par)
    }
  }
}
