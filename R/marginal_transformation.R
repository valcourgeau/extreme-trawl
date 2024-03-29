#' Transforms data from `GPD(1, 1 + kappa)` to `GPD(xi, sigma)`.
#' @param x Data to transform.
#' @param params_std Vector `c(xi, sigma, kappa)`.
#' @return Vector of `GPD(xi, sigma)` data.
#' @examples
#' xi <- 1
#' sig <- 3
#' kap <- 2
#' tmp <- evir::rgpd(100, xi = 1, beta = 1 + kap)
#' transformation_map(tmp, c(xi, sig, kap))
#' @export
transformation_map <- function(x, params_std) {
  # From GPD(1, 1 + kappa) to GPD(xi, sigma)

  params_trf <- parametrisation_translator(
    params = params_std, parametrisation = "standard", target = "transform"
  )
  # from GPD(xi, sigma) to Unif(0, 1)
  integral_trf <- evir::pgpd(
    q = x[which(x > 0.0)], mu = 0.0, xi = params_trf[1], beta = params_trf[2]
  )
  x[which(x > 0.0)] <- evir::qgpd(
    p = integral_trf, mu = 0.0, xi = params_std[1], beta = params_std[2]
  )

  return(x)
}

#' Inverse transformation function from `GPD(xi, sigma)` to `GPD(1, 1 + kappa)`.
#' @param x Data to transform
#' @param params_std Standard parametrisation parameters `(xi, sigma)`.
#' @return `GPD(1, 1 + kappa)` data.
#' @examples
#' n <- 10
#' xi <- .1
#' sig <- .5
#' kap <- 19.
#' data <- evir::rgpd(n, xi = xi, beta = sig)
#' transformation_map_inverse(data, c(xi, sig, kap))
#' @export
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


#' Jacobian of the transformation function.
#' @param params_std Standard parametrisation parameters `(xi, sigma)`.
#' @return Jacobian of the transformation function.
#' @examples
#' xi <- .1
#' sig <- .5
#' kap <- 19.
#' transformation_jacobian(c(xi, sig, kap))
#' @export
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
