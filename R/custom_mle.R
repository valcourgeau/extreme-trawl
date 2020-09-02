GetInitialGuessAndBounds <- function(data, max_length=20000, minus_mult=0.5, plus_multiplier=1.5){
  length_data <- min(max_length, length(data))
  data <- data[1:length_data]

  init_guess <- as.numeric(evir::gpd(data, threshold = 0, method="pwm")$par.ests) # (xi, sigma)

  p <- length(which(data>0)) / length(data)
  kap <- (1 - p) / p
  init_guess <- c(init_guess, kap)
  minus_init <- minus_mult * init_guess
  plus_init <- plus_multiplier * init_guess
  bds <- cbind(minus_init, plus_init)
  lower <- apply(bds, 1, min)
  upper <- apply(bds, 1, max)

  return(list(init_guess=init_guess, lower=lower, upper=upper))
}

CustomLikelihood <- function(data){
  # with 0 and non-0
  # returns a function of parameters

  return(function(par){
    data_for_mle <- data
    p <- length(which(data_for_mle>0)) / length(data_for_mle)

    kap <- par[3]
    p_zero <- 1 - 1. / (1. + kap)

    like <- evir::dgpd(data_for_mle[data_for_mle>0.0], beta = par[2], xi=par[1], mu=0.0) * (1-p_zero)
    like <- like[!is.na(like)] + 1e-7
    log.like <- sum(log(like)) + length(which(data_for_mle == 0.0)) * log(p_zero)

    log.like <- min(log.like, 1e9)
    log.like <- max(log.like, -1e9)

    return(-log.like)
  })
}


CustomMarginalMLE <- function(data, parametrisation='standard'){
  get_init_guess <- GetInitialGuessAndBounds(data)
  init_guess <- get_init_guess$init_guess
  lower <- get_init_guess$lower
  upper <- get_init_guess$upper
  fn_mle <- CustomLikelihood(data = data)

  res <- stats::optim(par = init_guess[1:3], fn_mle, method='L-BFGS-B', lower=lower, upper=upper)$par
  return(res)
}


CompositeLikelihood <- function(data){
  # implements custom likelihood with jacobian
  # returns a function of parameters
  return(function(par){
    data_for_mle <- data

    p_non_zero <- 1. / (1. + par[3])

    trf_data <- TransformationMapInverse(data_for_mle, params_std=par)
    custom_lk <- CustomLikelihood(data = trf_data) # needs to inject transformed data here
    trf_lk <- custom_lk(par = c(1., 1. + par[3], par[3])) # GPD(1, 1+kappa)

    jacob <- TransformationJacobian(params_std = par)
    jacob_vals <- jacob(data_for_mle)
    jacob_vals <- -sum(log(jacob_vals[!is.na(jacob_vals)] + 1e-7))

    composite_lk <- jacob_vals + trf_lk
    return(composite_lk)
  })
}


CompositeMarginalMLE <- function(data){
  get_init_guess <- GetInitialGuessAndBounds(data)
  init_guess <- get_init_guess$init_guess
  lower <- get_init_guess$lower
  upper <- get_init_guess$upper
  fn_mle <- CompositeLikelihood(data = data)

  return(stats::optim(par = init_guess[1:3], fn_mle, method='L-BFGS-B', lower=lower, upper=upper, control = list())$par)
}


CompositeLikelihoodScore <- function(params, max_length=100){
  return(
    function(data){
      n_row <- min(max_length, length(data))
      data <- data[1:n_row]
      grad_composite <- t(vapply(data, function(x){
        composite_lk <- CompositeLikelihood(x)
        return(pracma::grad(composite_lk, x0 = params))
    }, rep(0, length(params))))
      print(dim(grad_composite))
      return(grad_composite)
    }) # data_length x length(params)
}

CompositeMarginalHAC <- function(data, params, k=10, max_length=100){
  lk_score <- CompositeLikelihoodScore(params, max_length)
  score_acf <- AutocovarianceMatrix(lk_score(data), params, k)
  return(MakeHAC(score_acf))
}

