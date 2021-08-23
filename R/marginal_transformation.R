#' Transforms data from `GPD(1, 1 + kappa)` to `GPD(xi, sigma)`.
#' @param x Data to transform.
#' @param params_std Vector `c(xi, sigma, kappa)`.
#' @return Vector of `GPD(xi, sigma)` data.
#' @examples
#' @export
transformation_map <- function(x, params_std) {
  # From GPD(1, 1 + kappa) to GPD(xi, sigma)

  params_trf <- parametrisation_translator(
    params = params_std,
    parametrisation = "standard",
    target = "transform"
  )
  # from GPD(xi, sigma) to Unif(0, 1)
  integral_trf <- evir::pgpd(
    q = x[which(x > 0.0)],
    mu = 0.0,
    xi = params_trf[1],
    beta = params_trf[2]
  )
  x[which(x > 0.0)] <- evir::qgpd(
    p = integral_trf,
    mu = 0.0,
    xi = params_std[1],
    beta = params_std[2]
  )

  return(x)
}

transformation_map_inverse <- function(x, params_std) {
  # From GPD(xi, sigma) to GPD(1, 1 + kappa)

  params_trf <- parametrisation_translator(
    params = params_std,
    parametrisation = "standard",
    target = "transform"
  )

  integral_trf <- evir::pgpd(
    q = x[which(x > 0.0)], mu = 0.0, xi = params_std[1], beta = params_std[2]
  )
  x[which(x > 0.0)] <- evir::qgpd(
    p = integral_trf, mu = 0.0, xi = params_trf[1], beta = params_trf[2]
  )

  return(x)
}

transformation_jacobian <- function(params_std) {
  # parametrisation should be standard

  params_trf <- parametrisation_translator(
    params = params_std,
    parametrisation = "standard",
    target = "transform"
  )

  jacob_function <- function(z) {
    # takes transformed inputs
    z <- as.vector(z)
    z_non_zero <- z[z > 0.0]
    upper <- evir::dgpd(
      x = z_non_zero, mu = 0, beta = params_std[2], xi = params_std[1]
    )
    inv_z <- transformation_map_inverse(x = z_non_zero, params_std = params_std)
    lower <- evir::dgpd(inv_z, mu = 0, beta = params_trf[2], xi = params_trf[1])

    final_val <- rep(1.0, length(z))
    final_val[z > 0.0] <- pmin(pmax(upper / lower, 1e-9), 1e9)

    return(final_val)
  }

  return(jacob_function)
}
