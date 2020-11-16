library('zeallot')

GridFoundations <- function(n, vanishing_depth, x=1){
  # returns a triangular matrix

  # n <- n + vanishing_depth# accomodate for warm-up

  n_elems_per_line <- vanishing_depth
  index_i <- unlist(lapply(1:max(2, n_elems_per_line), FUN = function(i){return(1:(n))}))
  index_j <- unlist(lapply(1:max(2,n_elems_per_line), FUN = function(i){return(i:(n+i-1))}))
  return(Matrix::sparseMatrix(i=index_i, j=index_j, use.last.ij = T, x = x))
}

# n <- 14
# vanishing_depth <- 4
# gf <- GridFoundations(n = 14, vanishing_depth = vanishing_depth)
# gf[1:14,]
# k <- 13
#
# gf[max(c(1, k-vanishing_depth+1)):min(c(n,k)), k:min(c(n,vanishing_depth+k-1))]

TrawlSlicing <- function(n, vanishing_depth, trawl_parameter, type='exp'){
  # returns matrix with
  # [[S(1,1), S(2,1), S(3,1), ..., S(vanishing_depth,1)],
  #  [S(2,2), S(3,2), S(4,2), ..., S(vanishing_depth+1,2)],
  #  ...,
  #  [S(n,n), S(3,2), S(4,2), ..., S(vanishing_depth+n,n)],
  # ]

  b_funcs <- GetTrawlFunctions(type=type)
  B1_func <- b_funcs[[1]]
  B2_func <- b_funcs[[2]]
  B3_func <- b_funcs[[3]]

  one_split <- diff(-B2_func(param = trawl_parameter, h = 0:vanishing_depth))
  slices <- matrix(rep(one_split, n), ncol = vanishing_depth, byrow = T)
  area_A <- B2_func(trawl_parameter, 0.0)
  slices[2:nrow(slices),] <- slices[2:nrow(slices),] - rep(1, n-1) %o% c(one_split[2:length(one_split)], 0.0)

  return(slices/area_A) # divide by \mu^{leb}(A)
}

# tmp <- TrawlSlicing(n = 5, vanishing_depth = 10, trawl_parameter = 0.3)

GammaGrid <- function(alpha, beta, n, vanishing_depth, trawl_parameter, type='exp'){
  n_diags <- max(2, vanishing_depth)
  n_non_zero_elements <- n_diags*n

  gamma_shapes <- GenerateShapes(alpha = alpha,
                                 n = n,
                                 vanishing_depth = vanishing_depth,
                                 trawl_parameter = trawl_parameter,
                                 type = type)
  stopifnot(length(gamma_shapes) == n_non_zero_elements)
  gamma_sim_vals <- rgamma(n = n_non_zero_elements,
                           shape = gamma_shapes,
                           rate = beta)
  return(GridFoundations(n = n, vanishing_depth = vanishing_depth, x = gamma_sim_vals))
}

GenerateShapes <- function(alpha, n, vanishing_depth, trawl_parameter, type='exp'){
  stopifnot(length(alpha) == 1)

  slices <- TrawlSlicing(n = n, vanishing_depth = vanishing_depth, trawl_parameter = trawl_parameter, type = type) # already standardised
  return(alpha*as.vector(slices))
}

BlockIndex <- function(n_block, n, vanishing_depth){
  return(
    list(
      block_i=max(c(1, n_block-vanishing_depth+1)):min(c(n,n_block)),
      block_j=n_block:min(c(n,vanishing_depth+n_block-1)))
  )
}

GammaOrchestra <- function(scaled_gamma_grid, parallel=T){
  n <- dim(scaled_gamma_grid)[1]
  vanishing_depth <- dim(scaled_gamma_grid)[2] - n + 1

  if(parallel){
    cores <- detectCores(logical = TRUE)
    cl <- makeCluster(cores)
    clusterExport(cl, c('BlockIndex', 'n', 'vanishing_depth'))
    clusterEvalQ(cl, library(Matrix))
    tmp <- parallel::parLapply(cl = cl,
                               X = n:1,
                               fun = function(i){
                                 blck_ind <- BlockIndex(i, n = n, vanishing_depth = vanishing_depth)
                                 return(sum(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j]))
                               })
    parallel::stopCluster(cl)
    return(unlist(tmp))
  }else{
    return(vapply(1:n, FUN = function(i){
      blck_ind <- BlockIndex(i, n = n, vanishing_depth = vanishing_depth)
      return(sum(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j]))
    }, FUN.VALUE = 1.0))
  }
}

# n <- 10
# vanishing_depth <- 3
# gf <- GridFoundations(n, vanishing_depth)
#
# tmp <- GammaGrid(alpha = 3, beta = 20, n = 2500, vanishing_depth = 30, trawl_parameter = 0.1)
# tmp
# plot(GammaOrchestra(tmp))
# #
# b_funcs <- GetTrawlFunctions(type=type)
# B1_func <- b_funcs[[1]]
# B2_func <- b_funcs[[2]]
# B3_func <- b_funcs[[3]]
# B2_func(0.3, 0:10)
#
# (-B2_func(param = 0.8, h = 0:k)) %>% diff

PrintVanishingCoverage <- function(trawl_parameter, vanishing_depth, type='exp', get_value=F){
  b_funcs <- GetTrawlFunctions(type=type)
  B1_func <- b_funcs[[1]]
  B2_func <- b_funcs[[2]]
  B3_func <- b_funcs[[3]]

  max_val <- B2_func(trawl_parameter, h = 0.0)

  coverage <- -sum(diff(abs((B2_func(param = trawl_parameter, h = c(0,vanishing_depth))))))/abs(max_val)
  if(get_value){
    return(coverage)
  }else{
    cat('Coverage:', round(coverage*100.0, 2), '%\n')
  }
}

# PrintVanishingCoverage(trawl_parameter = 0.3, vanishing_depth = 10)

TrawlSimulation <- function(alpha, beta, n, vanishing_depth, trawl_parameter, type, parallel=F){
  gamma_grid <- GammaGrid(alpha = alpha,
                          beta = beta,
                          n = n,
                          vanishing_depth = vanishing_depth,
                          trawl_parameter = trawl_parameter,
                          type = type)
  return(GammaOrchestra(gamma_grid, parallel = parallel))
}

# profvis::profvis({
#   set.seed(42)
#   trawl_sim <- TrawlSimulation(3,3,n=5000,vanishing_depth=30,0.3,'exp',parallel = F)
# })
# profvis::profvis({
#   set.seed(42)
#   trawl_sim <- TrawlSimulation(3,3,n=5000,vanishing_depth=30,0.3,'exp',parallel = T)
# })
#
# trawl_sim <- TrawlSimulation(3,3,n=10000,vanishing_depth=30,0.3,'exp')
# trawl_sim %>% plot
# trawl_sim %>% acf
# lines(1:15, exp(-0.3*1:15))
#
# plot(density(trawl_sim))
# lines(0:50/10, dgamma(0:50/10, shape = 3, rate = 3))

ExceedancesSimulation <- function(params, n, vanishing_depth, type, m=1, parametrisation='standard', parallel=F, algo='standard'){
  #TODO have a function that does exactly that: given a parametrisation type, changes it to noven
  # set.seed(42)
  params_trawl <- ParametrisationTranslator(params = params[1:3], parametrisation = parametrisation, target = 'transform')
  kappa <- params[3]
  trawl_parameter <- params[4:length(params)]

  # print(
  #   fitdistrplus::fitdist(trawl_simulation, distr = "gamma", method = "mle")
  # )

  coverage <- PrintVanishingCoverage(trawl_parameter = trawl_parameter, vanishing_depth = vanishing_depth, type=type, get_value = T)

  if(parallel){
    cores <- parallel::detectCores(logical = TRUE)
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(
      cl, c('TransformationMapInverse',
            'TransformationMap',
            'TransformationJacobian',
            'ParametrisationTranslator',
            'PairwiseLikelihood',
            'CompositeMarginalMLE',
            'TrawlGMM',
            'TrawlAutocorrelation',
            'EVTrawlFit',
            'TrawlSimulation',
            'PrintVanishingCoverage',
            'GammaOrchestra',
            'BlockIndex',
            'GridFoundations',
            'TrawlSlicing',
            'GammaGrid',
            'GenerateShapes',
            'ExceedancesSimulation',
            GetTrawlEnvsList()))

    sim_fn <- function(i){
      ExceedancesSimulation(
        params = params, n = n,
        vanishing_depth = vanishing_depth, type = type,
        m = 1, parametrisation = parametrisation,
        parallel = F, algo = algo)
    }

    parallel::clusterExport(
      cl, c('n', 'm', 'vanishing_depth', 'type'))


    sims <- parallel::parLapply(
      cl = cl,
      X = 1:m,
      fun = sim_fn
    )
    parallel::stopCluster(cl)
    return(sims)
  }


  # corr_uniform <- pgamma(trawl_simulation, shape=alpha, rate=beta)
  # TODO should be 1-exp()

  if(algo == 'standard'){
    print('standard algo')
    trawl_simulation <- TrawlSimulation(alpha = 1.,
                                        beta = 1.,
                                        trawl_parameter = trawl_parameter,
                                        n = n,
                                        vanishing_depth = vanishing_depth,
                                        type = type,
                                        parallel = parallel)
    probabilities_zero <- 1-exp(-kappa*trawl_simulation)
    uniform_samples <- runif(n = n, min = 0, max = 1)
    who_is_extreme <- uniform_samples > probabilities_zero

    exceedances <- rep(0, n)
    exceedances[who_is_extreme] <- rexp(n = length(which(who_is_extreme)), rate = trawl_simulation[who_is_extreme])
    exceedances[who_is_extreme] <- TransformationMap(
      x = exceedances[who_is_extreme],
      params_std = params[1:3]
    )

    return(list(exceedances=exceedances, latent=trawl_simulation, coverage=coverage))
  }else{
    if(algo == 'cross'){
      print('cross algo')


      sim_data <- lapply(
        1:m,
        function(i){
          ExceedancesSimulation(params = params, parametrisation = parametrisation, n = n, m = m, vanishing_depth = vanishing_depth, type = type, parallel = F, algo = 'standard')
        }
      )

      backbone_data <- ExceedancesSimulation(params = params, parametrisation = parametrisation, n = n, m = m, vanishing_depth = vanishing_depth, type = type, parallel = F, algo = 'standard')
      backbone_exc <- backbone_data$exceedances

      sim_data <- abind::abind(lapply(sim_data, function(x){do.call(cbind, x)}), along=3) # [,2,] for latent; [,1,] for exceedances
      sim_exc <- sim_data[,1,]
      sim_latent <- sim_data[,2,]

      exceedances <- rep(0, n)
      print(apply(sim_exc[which(backbone_exc > 0),], 1, function(x){mean(x[x>0])}))
      exceedances[which(backbone_exc > 0)] <- unlist(apply(sim_exc[which(backbone_exc > 0),], 1, function(x){mean(x[x>0])}))

      acf(exceedances)
      return(list(exceedances=exceedances, coverage=coverage))
    }
    if(algo == 'corr_unif'){
      trawl_simulation <- TrawlSimulation(alpha = 1.,
                                          beta = 1.,
                                          trawl_parameter = trawl_parameter,
                                          n = n,
                                          vanishing_depth = vanishing_depth,
                                          type = type,
                                          parallel = parallel)
      trawl_simulation_unif <- TrawlSimulation(alpha = 1.,
                                          beta = 1.,
                                          trawl_parameter = trawl_parameter,
                                          n = n,
                                          vanishing_depth = vanishing_depth,
                                          type = type,
                                          parallel = parallel)

      probabilities_zero <- 1-exp(-kappa*trawl_simulation)
      corr_uniform <- pgamma(trawl_simulation_unif, shape=1, rate=1)

      who_is_extreme <- corr_uniform > probabilities_zero

      exceedances <- rep(0, n)
      exceedances[who_is_extreme] <- rexp(n = length(which(who_is_extreme)), rate = trawl_simulation[who_is_extreme])
      exceedances[who_is_extreme] <- TransformationMap(
        x = exceedances[who_is_extreme],
        params_std = params[1:3]
      )

      return(list(exceedances=exceedances, latent=trawl_simulation, coverage=coverage))
    }
    if(algo == 'dynamic_latent' | algo == 'dynamic_uniform'){
      trawl_simulation <- TrawlSimulation(alpha = 1.,
                                          beta = 1.,
                                          trawl_parameter = trawl_parameter,
                                          n = n,
                                          vanishing_depth = vanishing_depth,
                                          type = type,
                                          parallel = parallel)

      b_funcs <- GetTrawlFunctions(type=type)
      B1_func <- b_funcs[[1]]
      B2_func <- b_funcs[[2]]
      B3_func <- b_funcs[[3]]
      params_dynamic_thres <- c(1, 1, kappa, trawl_parameter)

      A_value <- B1_func(param=trawl_parameter, h=1) + B2_func(param=trawl_parameter, h=1)
      b_1_minus_0 <- - B1_func(param=trawl_parameter, h=1) / A_value
      b_0_minus_1 <- - B3_func(param=trawl_parameter, h=1) / A_value
      b_0_1 <- - B2_func(param=trawl_parameter, h=1) / A_value

      probabilities_zero <- 1-exp(-kappa*trawl_simulation)
      # corr_uniform <- pgamma(trawl_simulation_unif, shape=1, rate=1)
      uniform_samples <- runif(n = n, min = 0, max = 1)

      prev_sample <- NULL
      exceedances <- apply(
        cbind(probabilities_zero, uniform_samples, trawl_simulation),
        MARGIN = 1,
        FUN = function(p_u_t){
          if(is.null(prev_sample)){
            extreme_proba <- 1/(1+kappa)
          }else{
            extreme_proba <- (1+kappa+prev_sample)^{1+b_0_minus_1} *
              (1+2*kappa+prev_sample)^{b_0_1} *
              (1+kappa)^{b_1_minus_0}
            if(algo == 'dynamic_latent'){
              if(prev_sample == 0.0){
                extreme_proba <-  (1/(1+kappa) - extreme_proba/(1+kappa)) / (1-1/(1+kappa)) # /(1+kappa)
              }
            }else{
              if(algo == 'dynamic_uniform'){
                if(prev_sample == 0.0){
                  extreme_proba <- (1/(1+kappa) - extreme_proba/(1+kappa)) / (1-1/(1+kappa))
                }
              }
            }
          }
          if(algo == 'dynamic_latent'){
            if(1-p_u_t[1] > extreme_proba){
              prev_sample <<- rexp(n=1, rate=p_u_t[3])
            }else{
              prev_sample <<- 0.0
              # prev_sample <<- rexp(n=1, rate=p_u_t[3])
            }
          }else{
            if(algo == 'dynamic_uniform'){
              if(p_u_t[2] > extreme_proba){
                prev_sample <<- 0.0
              }else{
                prev_sample <<- rexp(n=1, rate=p_u_t[3])
              }
            }
          }


          return(prev_sample)
        }
      )

      who_is_extreme <- exceedances > 0

      # exceedances <- rep(0, n)
      # exceedances[who_is_extreme] <- rexp(n = length(which(who_is_extreme)), rate = trawl_simulation[who_is_extreme])
      exceedances[who_is_extreme] <- TransformationMap(
        x = exceedances[who_is_extreme],
        params_std = params[1:3]
      )

      return(list(exceedances=exceedances, latent=trawl_simulation, coverage=coverage))
    }
  }

  stop('Wrong simulation algo.')
}

