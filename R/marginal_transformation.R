TransformationMap <- function(x, params_std){
  # From GPD(1, 1 + kappa) to GPD(xi, sigma)

  params_trf <- ParametrisationTranslator(
      params = params_std,
      parametrisation = "standard",
      target = "transform")
  # from GPD(xi, sigma) to Unif(0, 1)
  integral_trf <- evir::pgpd(q = x[which(x>0.0)],
                       mu = 0.0,
                       xi = params_std[1],
                       beta = params_std[2])
  x[which(x>0.0)] <- evir::qgpd(p = integral_trf,
                                mu = 0.0,
                                xi = params_trf[1],
                                beta = params_trf[2])
  return(x)
}

TransformationMapInverse <- function(x, params_std){
  # From GPD(xi, sigma) to GPD(1, 1 + kappa)

  params_trf <- ParametrisationTranslator(
      params = params_std,
      parametrisation = "standard",
      target = "transform")

  integral_trf <- evir::pgpd(q = x[which(x>0.0)],
                             mu = 0.0,
                             xi = params_trf[1],
                             beta = params_trf[2])

  x[which(x>0.0)] <- evir::qgpd(p = integral_trf,
                                mu = 0.0,
                                xi = params_std[1],
                                beta = params_std[2])
  return(x)
}

TransformationJacobian <- function(params_std){
  #parametrisation should be standard

  params_trf <- ParametrisationTranslator(
    params = params_std,
    parametrisation = "standard",
    target = "transform")

  jacob_function <- function(z){
    #takes transformed inputs
    z <- as.vector(z)
    z_non_zero <- z[z>0.0]
    upper <- evir::dgpd(z_non_zero, mu = 0, beta = params_std[2], xi = params_std[1])
    inv_z <- TransformationMapInverse(z_non_zero, params_std = params_std)
    lower <- evir::dgpd(inv_z, mu = 0, beta = params_trf[2], xi = params_trf[1])

    final_val <- rep(1.0, length(z))
    final_val[z>0.0] <- pmin(pmax(upper/lower, 1e-9), 1e9)

    return(final_val)
  }

  return(jacob_function)
}

NovenJacobian <- function(params_std){
  stop('Not Implemented')
}


#
#
# set.seed(42)
# gpd_sample_1 <- evir::rgpd(n = 10000, mu = 0.0, beta = (1+1)/1., xi = 1.)
# plot(log(sort(gpd_sample_1)), log(sort(TransformationJacobian(c(1.0, 1, 1))(gpd_sample_1))))
#
# set.seed(42)
# gpd_sample <- evir::rgpd(n = 10000, mu = 0.0, beta = 5, xi = -0.2)
# tmp <- TransformationMap(
#     TransformationMapInverse(gpd_sample, params = c(-0.2, 5, 3)),
#     params = c(-0.2, 5, 3))
# plot(abs(sort(tmp)-sort(gpd_sample)) < 1e-10, main='one is good')
# line(0:20, 0:20)
# evir::gpdFit(tmp, 0.0)
