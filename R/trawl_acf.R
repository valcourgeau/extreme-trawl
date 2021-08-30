.onUnload <- function(libpath) { # nolint
  library.dynam.unload("extreme.trawl", libpath)
}

#' R-only Single Trawl ACF calculations.
#' @param h Single ACF horizon.
#' @param alpha Gamma `alpha` (or shape) parameter.
#' @param betal Gamma `beta` (or rate) parameter.
#' @param rho Trawl parameter(s).
#' @param kappa Extreme frequency control parameter.
#' @param delta Discrete time step.
#' @param end_seq Horizon cut-off for integral approximation.
#' @param type Trawl type (e.g. `"exp`, "`sum_exp`", etc.).
#' @param cov Boolean, correlation or covariance (if `FALSE`).
#'     Defaults to correlation.
#' @return Trawl autocorrelation values with R-only functions.
#' @export
acf_trawl_single <- function(h, alpha, beta, rho, kappa,
                             delta = 0.1, end_seq = 50,
                             type = "exp", cov = F) {
  # Compute ACF with trawl process as latent
  seq_kappa <- seq(kappa, kappa + end_seq, by = delta)
  trawl_fct <- get_trawl_functions(type)
  b_1_func <- trawl_fct[[1]]
  b_2_func <- trawl_fct[[2]]
  b_3_func <- trawl_fct[[3]]

  a_total <- b_1_func(param = rho, h = h) + b_2_func(param = rho, h = h)
  b_h_minus_0 <- -alpha * b_1_func(param = rho, h = h) / a_total
  b_0_minus_h <- -alpha * b_3_func(param = rho, h = h) / a_total
  b_0_h <- -alpha * b_2_func(param = rho, h = h) / a_total

  sum_over_x <- vapply(seq_kappa,
    FUN = function(x) {
      x <- x + delta / 2
      sum_over_y <- vapply(
        X = seq_kappa, FUN = function(y) {
          y <- y + delta / 2
          tmp1 <- (1 + x / beta)^b_h_minus_0
          tmp2 <- (1 + (x + y) / beta)^b_0_h
          tmp3 <- (1 + y / beta)^b_0_minus_h
          tmp4 <- (1 + (x + y) / beta)^{
            -alpha
          }
          return(c(tmp1 * tmp2 * tmp3, tmp4))
        },
        FUN.VALUE = rep(0, 2)
      )
      total_sum_over_y <- apply(sum_over_y, MARGIN = 1, sum)
      total_x <- (1 + x / beta)^{
        -alpha
      }
      return(c(total_sum_over_y, total_x))
    },
    FUN.VALUE = rep(0, 3)
  )

  final_sum <- apply(sum_over_x, MARGIN = 1, sum)
  final_sum <- final_sum * delta
  final_sum[1:2] <- final_sum[1:2] * delta

  res <- final_sum[1]
  res_0 <- final_sum[2]
  first_mom_sq <- final_sum[3]^2

  if (cov) {
    return(res - first_mom_sq)
  } else {
    return((res - first_mom_sq) / (res_0 - first_mom_sq))
  }
}

#' Rcpp-accelerated Single Trawl ACF calculations.
#' @param h ACF horizon.
#' @param alpha Gamma `alpha` (or shape) parameter.
#' @param betal Gamma `beta` (or rate) parameter.
#' @param rho Trawl parameter(s).
#' @param kappa Extreme frequency control parameter.
#' @param delta Discrete time step.
#' @param end_seq Horizon cut-off for integral approximation.
#' @param type Trawl type (e.g. `"exp`, "`sum_exp`", etc.).
#' @param cov Boolean, correlation or covariance (if `FALSE`).
#'     Defaults to correlation.
#' @return Trawl autocorrelation values with Cpp functions.
#' @export
cpp_acf_trawl_single <- function(h, alpha, beta, rho,
                                 kappa, delta = 0.1,
                                 end_seq = 50, type = "exp", cov = F) {
  # Compute ACF with trawl process as latent
  seq_kappa <- seq(kappa, kappa + end_seq, by = delta)
  trawl_fct <- get_trawl_functions(type)
  b_1_func <- trawl_fct[[1]]
  b_2_func <- trawl_fct[[2]]
  b_3_func <- trawl_fct[[3]]

  a_total <- b_1_func(param = rho, h = h) + b_2_func(param = rho, h = h)
  b_0_minus_h <- -alpha * b_3_func(param = rho, h = h) / a_total
  b_0_h <- -alpha * b_2_func(param = rho, h = h) / a_total

  res_0 <- square_moment(
    xs = seq_kappa, delta = delta, beta = beta,
    b_oh = b_0_h, b_o_exc_h = b_0_minus_h
  )
  res <- cross_moment(
    xs = seq_kappa, delta = delta, beta = beta,
    b_oh = b_0_h, b_o_exc_h = b_0_minus_h
  )
  first_mom <- first_moment(
    xs = seq_kappa, delta = delta, beta = beta,
    b_oh = b_0_h, b_o_exc_h = b_0_minus_h
  )
  first_mom_sq <- first_mom^2

  if (cov) {
    return(res - first_mom_sq)
  } else {
    return((res - first_mom_sq) / (res_0 - first_mom_sq))
  }
}

cross_moment_trawls <- function(h, alpha, beta, rho,
                                kappa, delta = 0.1,
                                end_seq = 50, type = "exp") {
  seq_kappa <- seq(kappa, kappa + end_seq, by = delta)
  trawl_fct <- get_trawl_functions(type)
  b_1_func <- trawl_fct[[1]]
  b_2_func <- trawl_fct[[2]]
  b_3_func <- trawl_fct[[3]]

  a_total <- b_1_func(param = rho, h = h) + b_2_func(param = rho, h = h)
  b_h_minus_0 <- -alpha * b_1_func(param = rho, h = h) / a_total
  b_0_minus_h <- -alpha * b_3_func(param = rho, h = h) / a_total
  b_0_h <- -alpha * b_2_func(param = rho, h = h) / a_total

  sum_over_x <- vapply(seq_kappa,
    FUN = function(x) {
      x <- x + delta / 2
      sum_over_y <- vapply(
        X = seq_kappa, FUN = function(y) {
          y <- y + delta / 2
          tmp1 <- (1 + x / beta)^b_h_minus_0
          tmp2 <- (1 + (x + y) / beta)^b_0_h
          tmp3 <- (1 + y / beta)^b_0_minus_h
          return(tmp1 * tmp2 * tmp3)
        },
        FUN.VALUE = rep(0, 1)
      )
      return(sum(sum_over_y))
    },
    FUN.VALUE = rep(0, 1)
  )

  final_sum <- sum(sum_over_x)
  res <- final_sum * delta^2
  return(res)
}

#' R-only Trawl ACF function approximation.
#' @param h ACF horizon.
#' @param alpha Gamma `alpha` (or shape) parameter.
#' @param betal Gamma `beta` (or rate) parameter.
#' @param rho Trawl parameter(s).
#' @param kappa Extreme frequency control parameter.
#' @param delta Discrete time step.
#' @param end_seq Horizon cut-off for integral approximation.
#' @param type Trawl type (e.g. `"exp`, "`sum_exp`", etc.).
#' @param cov Boolean, covariance (if `TRUE`) or correlation.
#'     Defaults to covariance.
#' @return Trawl autocorrelation values with R-only functions.
#' @export
acf_trawl <- function(h, alpha, beta, kappa,
                      rho, delta = 0.5,
                      type = "exp", cov = T) {
  trawl_acf_fn <- function(h) {
    acf_trawl_single(
      h,
      alpha = alpha, beta = beta, kappa = kappa,
      rho = rho, delta = delta, type = type, cov = cov
    )
  }
  return(vapply(h, trawl_acf_fn, 1))
}

acf_trawl_revisited_num_approx <- function(h,
                                           alpha, beta,
                                           kappa, rho,
                                           delta = 0.5,
                                           type = "exp",
                                           cov = T) {
  trawl_acf_fn <- function(h) {
    acf_trawl_single(
      h,
      alpha = alpha, beta = beta, kappa = kappa,
      rho = rho, delta = delta, type = type, cov = cov
    )
  }
  cores <- parallel::detectCores(logical = TRUE)
  cl <- parallel::makeCluster(cores - 1)
  parallel::clusterCall(
    cl, c(
      "acf_trawl", "square_moment", "cross_moment", "first_moment",
      "trawl_acf_fn", get_trawl_envs_list()
    )
  )
  acf_vals <- parallel::parLapply(cl, X = h, fun = trawl_acf_fn)
  parallel::stopCluster(cl)
  return(acf_vals)
}

#' Rcpp-accelerated Trawl ACF function approximation.
#' @param h ACF horizon.
#' @param alpha Gamma `alpha` (or shape) parameter.
#' @param betal Gamma `beta` (or rate) parameter.
#' @param rho Trawl parameter(s).
#' @param kappa Extreme frequency control parameter.
#' @param delta Discrete time step.
#' @param end_seq Horizon cut-off for integral approximation.
#' @param type Trawl type (e.g. `"exp`, "`sum_exp`", etc.).
#' @param cov Boolean, correlation or covariance (if `FALSE`).
#'     Defaults to correlation.
#' @return Trawl autocorrelation values with Cpp-accelerated functions.
#' @export
cpp_acf_trawl <- function(h, alpha, beta, kappa,
                          rho, delta = 0.5,
                          type = "exp", end_seq = 50,
                          cov = T) {
  trawl_acf_fn <- function(h) {
    cpp_acf_trawl_single(
      h,
      alpha = alpha, beta = beta, kappa = kappa, end_seq = end_seq,
      rho = rho, delta = delta, type = type, cov = cov
    )
  }
  return(vapply(X = h, FUN.VALUE = 1.0, FUN = trawl_acf_fn))
}
