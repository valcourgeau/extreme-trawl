TrawlObjective <- function(data, depth, parametrisation='standard'){
  function(pars){
    kappa <- GetKappa(data = data, params = pars, parametrisation = 'standard')
    pars <- c(pars[1:2], kappa, pars[3:length(pars)])
    noven_pars <- ParametrisationTranslator(params = pars, parametrisation = parametrisation, target = 'noven')
    return(function(trawl_params){
      acf_vals <- vapply(c(0.01, 1:(depth)), function(h){
        acf_trawl(h, alpha = noven_pars[1], beta = noven_pars[2], kappa = noven_pars[3],
                  rho = trawl_params, delta = 0.5, end_seq = 50)}, 1)
      sample_cross_mom <- acf(data, plot = F, lag.max = depth)$acf

      return(sum((acf_vals-sample_cross_mom)))
    }
    )
  }
}

GMMObjective <- function(data, depth, omega='id', parametrisation='standard'){
  composite <- CustomLikelihood(data = data,
                                parametrisation=parametrisation)
  trawl_objective <- TrawlObjective(data = data,
                                    depth = depth,
                                    parametrisation = parametrisation)

  return(function(par){
    grad_vec <- c(
      pracma::grad(composite, x0 = par[1:2]),
      pracma::grad(trawl_objective(par), x0 = par[3:length(par)])
    )

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


pollution_data <- read.csv('data/clean_pollution_data.csv')
to <- TrawlObjective(data = pollution_data[,2],
                     depth = 10,
                     parametrisation = 'standard')
max_depth <- 10000
custom_mle <- CustomMarginalMLE(pollution_data[1:max_depth,2])
kappa <- GetKappa(pollution_data[,2] ,params = custom_mle, parametrisation = 'standard')
custom_mle_kappa <- c(custom_mle, kappa)
c_mle_kappa_rho <- c(custom_mle_kappa, 0.2)
gmm_obj <- GMMObjective(data = pollution_data[1:max_depth,2], depth = 10)
gmm_obj(c_mle_kappa_rho[-3])


profvis::profvis({
  gmm_obj(custom_mle_kappa)
})

optim(gmm_obj, par = c(custom_mle, 0.2), method = 'L-BFGS-B',
      lower=c(1e-2, 1e-2, 1e-2),
      upper=c(0.5,2,1.0),control = list(trace=5))

profvis::profvis()
