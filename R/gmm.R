TrawlGMM <- new.env()

TrawlGMM$Diff <- function(x, y){
  return(sum(x - y))
}

TrawlGMM$TrawlObjective <- function(data, depth, type='exp', metric=measures::SSE){
  function(pars){
    # pars should be (xi, sigma, kappa, rho)
    noven_pars <- ParametrisationTranslator(params = pars[1:3], parametrisation = 'standard', target = 'transform')
    sample_cross_mom <- acf(data, plot = F, lag.max = depth)$acf

    return(function(trawl_params){
        acf_vals <- TrawlAutocorrelation$AcfTrawlCollection(
            h=c(0.01, 1:(depth)), alpha = noven_pars[1],
            beta = noven_pars[2], kappa = noven_pars[3],
            rho = trawl_params, delta = 0.1, cov = F, type=type)

        return(metric(acf_vals, sample_cross_mom))
      }
    )
  }
}

TrawlGMM$FullGMMObjective <- function(data, depth, omega='id', type='exp'){
  composite <- CompositeLikelihood(data = data)
  trawl_objective <- TrawlGMM$TrawlObjective(data = data,
                                    depth = depth,
                                    type=type,
                                    metric=TrawlGMM$Diff) # without derivative on the last one
  return(function(par){
    grad_vec <- c(
      pracma::grad(composite, x0 = par[1:3]) / length(data),
      trawl_objective(par)(par[4:length(par)]) # pracma::grad(trawl_objective(par), x0 = par[4:length(par)])
    )

    # grad_vec[4] <- grad_vec[4] * sum(abs(grad_vec[1:3]))
    # grad_vec <- grad_vec / (sum(abs(grad_vec)) + 1)
    # grad_vec <- grad_vec / length(data)

    if(omega == 'id'){
      omega <- diag(rep(1, length(par)))
    }else{

      if(omega == 'centered'){
        omega <- diag(rep(1, length(par)))
      }
    }

    return(t(grad_vec) %*% omega %*% grad_vec)
  })
}

TrawlGMM$TwoStageGMMObjective <- function(data, depth, type='exp', metric=measures::SSE){
  pars <- CompositeMarginalMLE(data = data)
  trawl_obj <- TrawlGMM$TrawlObjective(data, depth, type=type, metric=metric) # return a function of the whole set of params
  return(trawl_obj(pars)) # function of trawl parameters
}

TrawlGMM$TrawlObjectiveDatapoint <- function(xs, k, mean_data, type='exp', metric=measures::SSE){
  # metric can be measures::SSE or measures::MAE
  # return the trawl objective for each k and data points

  return(
    function(pars){
      # pars should be (xi, sigma, kappa, rho)
      noven_pars <- ParametrisationTranslator(params = pars[1:3], parametrisation = 'standard', target = 'transform')

      return(function(trawl_params){
          acf_vals <- TrawlAutocorrelation$AcfTrawlCollection(
            h=k, alpha = noven_pars[1],
            beta = noven_pars[2], kappa = noven_pars[3],
            rho = trawl_params, delta = 0.1, cov = F, type=type)

          xs_stack <- prod(xs) - mean_data^2

          return(metric(acf_vals, xs_stack))
        }
      )
    }
  )
}

TrawlGMM$CompositeAndTrawlObjectiveDatapoint <- function(xs, k, mean_data, length_data, omega='id', type='exp'){
  composite <- CompositeLikelihood(data = xs)
  trawl_objective <- TrawlGMM$TrawlObjectiveDatapoint(
    xs=xs,
    k=k,
    mean_data=mean_data,
    type=type,
    metric=TrawlGMM$Diff)
  return(function(par){
    grad_vec <- c(
      -pracma::grad(composite, x0 = par[1:3]),
      trawl_objective(par)(par[4:length(par)])
    )

    return(grad_vec)
  })
}

TrawlGMM$WrapperTrawlObjectiveDatapoint <- function(xs, k, mean_data, length_data, omega='id', type='exp'){
  trawl_objective <- TrawlGMM$TrawlObjectiveDatapoint(
    xs=xs,
    k=k,
    mean_data=mean_data,
    type=type,
    metric=TrawlGMM$Diff)
  return(function(par){
    grad_vec <-  trawl_objective(par)(par[4:length(par)])
    return(grad_vec)
  })
}


TrawlGMM$TrawlGMMScore <- function(params, depth, type='exp', max_length=100){
  # Full Score function
  return(
    function(data){
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]
      mean_data <- mean(data)

      score_per_depth <- lapply(
        1:depth,
        function(k){
          xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
          t(apply(xs_stack, MARGIN = 1,
                  FUN = function(xs){
                      full_gmm <- TrawlGMM$CompositeAndTrawlObjectiveDatapoint(xs, k, mean_data, n_sample)
                      return(full_gmm(params))
                   }))})
      return(score_per_depth)
    }) # list of depth items data_length x length(params)
}

TrawlGMM$TrawlGMMScorePartial <- function(params, depth, type='exp', max_length=100){
  # Full Score function
  return(
    function(data){
      n_sample <- min(max_length, length(data))
      data <- data[1:n_sample]
      mean_data <- mean(data)

      score_per_depth <- lapply(
        1:depth,
        function(k){
          xs_stack <- cbind(data[1:(n_sample-k)], data[(k+1):(n_sample)])
          apply(xs_stack, MARGIN = 1,
                  FUN = function(xs){
                    full_gmm <- TrawlGMM$WrapperTrawlObjectiveDatapoint(xs, k, mean_data, n_sample)
                    return(full_gmm(params))
                  })})
      return(score_per_depth)
    }) # list of depth items data_length x length(params)
}

TrawlGMM$TrawlGMMHAC <- function(data, params, depth, k=10, type='exp', max_length=100){
  lk_score <- TrawlGMM$TrawlGMMScore(params, depth, type, max_length)
  pl_score_per_depth <- lk_score(data)

  score_acf_autocov_mat <- lapply(
    pl_score_per_depth,
    function(pl_score){AutocovarianceMatrix(pl_score, params, k)})
  pl_hac <- lapply(score_acf_autocov_mat, function(autocov_mat){MakeHAC(autocov_mat, near.pd=F)})
  return(as.matrix(Matrix::nearPD(Reduce(`+`, pl_hac)/depth)$mat)) # sum across clusters
}

TrawlGMM$TrawlGMMHACPartial <- function(data, params, depth, k=10, type='exp', max_length=100){
  # only the trawl parameters
  lk_score <- TrawlGMM$TrawlGMMScorePartial(params, depth, type, max_length)
  pl_score_per_depth <- lk_score(data)

  trawl_params <- params[4:length(params)]

  score_acf_autocov_mat <- lapply(
    pl_score_per_depth,
    function(pl_score){AutocovarianceMatrix(pl_score, trawl_params, k)})
  pl_hac <- lapply(score_acf_autocov_mat, function(autocov_mat){MakeHAC(autocov_mat)})
  return(Reduce(`+`, pl_hac)/depth) # sum across clusters
}




