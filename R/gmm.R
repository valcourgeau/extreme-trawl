TrawlObjective <- function(data, depth, type='exp'){
  function(pars){
    # pars should be (xi, sigma, kappa, rho)
    noven_pars <- ParametrisationTranslator(params = pars[1:3], parametrisation = 'standard', target = 'transform')
    return(function(trawl_params){
        acf_vals <- AcfTrawlCollection(
            h=c(0.01, 1:(depth)), alpha = noven_pars[1],
            beta = noven_pars[2], kappa = noven_pars[3],
            rho = trawl_params, delta = 0.1, cov = F, type=type)
        sample_cross_mom <- acf(data, plot = F, lag.max = depth)$acf

        return(sum((acf_vals-sample_cross_mom)^2))
      }
    )
  }
}

FullGMMObjective <- function(data, depth, omega='id', type='exp'){
  composite <- CompositeLikelihood(data = data)
  trawl_objective <- TrawlObjective(data = data,
                                    depth = depth,
                                    type=type)
  return(function(par){
    grad_vec <- c(
      pracma::grad(composite, x0 = par[1:3]),
      pracma::grad(trawl_objective(par), x0 = par[4:length(par)])
    )

    grad_vec[4] <- grad_vec[4] * sum(abs(grad_vec[1:3]))
    grad_vec <- grad_vec / sum(abs(grad_vec))

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

TwoStageGMMObjective <- function(data, depth, type='exp'){
  pars <- CompositeMarginalMLE(data = data)
  trawl_obj <- TrawlObjective(data, depth, type=type) # return a function of the whole set of params
  return(trawl_obj(pars)) # function of trawl parameters
}
