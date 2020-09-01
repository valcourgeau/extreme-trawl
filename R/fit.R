
EVTrawlFit <- function(data, depth, method, mode='two-stage', type='exp', bounds='config', cl=NULL, ...){
  stopifnot(mode %in% c('two-stage', 'full'))
  stopifnot(method %in% c('GMM', 'PL'))
  stopifnot(bounds %in% c('config', 'multiplier'))

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
        optim_fn <- TrawlGMM$TwoStageGMMObjective(data = data, depth = depth)
      }
    }

    # choose bounds
    lower <- trawl_cfg$lower
    upper <- trawl_cfg$upper
    if(bounds == 'multiplier'){
      init_trawl <- PairwiseLikelihood$InitGuess(data=data, depth=depth, n_trials=20, type=type)
      lower <- init_trawl * 0.5
      upper <- init_trawl * 1.5
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
          optim_fn <- TrawlGMM$FullGMMObjective(
            data = data, depth = depth)
        }
      }

      # choose bounds
      lower <- trawl_cfg$lower
      upper <- trawl_cfg$upper
      if(bounds == 'multiplier'){
        init_trawl <- PairwiseLikelihood$InitGuess(data=data, depth=depth, n_trials=20, type=type)
        lower <- init_trawl * 0.5
        upper <- init_trawl * 1.5
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

SubSampleFit <- function(data, sample.length, depth, method, mode='two-stage', type='exp', bounds='config',
                         trials, parallel=F, seed=42){
  # method 'PL' or 'GMM'
  # depth for PL is the length of blocks, GMM depth is the ACF depth
  n <- length(data)
  set.seed(seed)
  start_points <- sample(1:(n-sample.length), size = trials, replace = F)

  results_list <- list()
  results <- matrix(0, nrow=trials, ncol=3+GetTrawlParamsConfig(type)$n_params)

  if(parallel){
    cores <- parallel::detectCores(logical = TRUE)
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(
      cl, c('TransformationMapInverse',
            'TransformationMap',
            'TransformationJacobian',
            'ParametrisationTranslator',
            'PairwiseLikelihood',
            'TrawlGMM',
            'TrawlAutocorrelation',
            'EVTrawlFit',
            GetTrawlEnvsList()))

    sub_sample_time <- Sys.time()
    if(method == 'GMM'){
      # we perform parallel estimations
      results <- parallel::parLapply(
        X = start_points,
        cl = cl,
        fun = function(start_pt){
          EVTrawlFit(
            data = data[start_pt:(start_pt+sample.length)],
            depth = depth,
            method = method,
            type = type,
            bounds = bounds,
            cl = NULL)
         }
      )
    }else{
      # we perform PL computation in parallel
      results <- lapply(
        X = start_points,
        FUN = function(start_pt){
          EVTrawlFit(
            data = data[start_pt:(start_pt+sample.length)],
            depth = depth,
            method = method,
            type = type,
            bounds = bounds,
            cl = cl)
        }
      )
    }

    print(Sys.time() - sub_sample_time)
    parallel::stopCluster(cl)
    results <- matrix(unlist(results), ncol=length(results[[1]]), byrow = T)
  }else{
    print('No parallel trials but parallel PL.')
    sub_sample_time <- Sys.time()
    results<- t(vapply(start_points,
                       FUN = function(start_pt){
                         tmp_pr <- EVTrawlFit(data = data[start_pt:(start_pt+sample.length)],
                                              depth = depth,
                                              parametrisation = parametrisation,
                                              type = type,
                                              method = method,
                                              bounds = bounds,
                                              cl = NULL)
                       },
                       FUN.VALUE = rep(0, ncol(results))))
    print(Sys.time() - sub_sample_time)
  }

  results_list$estimators <- results
  results_list$sample.length <- sample.length
  results_list$depth <- depth
  results_list$mode <- mode
  results_list$method <- method
  results_list$bounds <- bounds
  results_list$type <- type
  results_list$data.length <- n
  results_list$trials <- trials
  results_list$start.indices <- start_points

  return(results_list)
}