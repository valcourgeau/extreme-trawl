.onUnload <- function (libpath) { library.dynam.unload("gammaextremes", libpath)}

acf_trawl <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp', cov=F){
  # Compute ACF with trawl process as latent
  seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
  trawl_fct <- GetTrawlFunctions(type)
  B1_func <- trawl_fct[[1]]
  B2_func <- trawl_fct[[2]]
  B3_func <- trawl_fct[[3]]

  b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))

  res <- 0
  first_mom <- 0
  res_0 <- 0
  sum_over_x <- vapply(seq_kappa, FUN = function(x){
    x <- x + delta / 2
    sum_over_y <- vapply(X = seq_kappa, FUN = function(y){
      y <- y + delta / 2
      tmp1 <- (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
      tmp2 <- (1+(x+y)/beta)^{-alpha}
      return(c(tmp1, tmp2))
    },
    FUN.VALUE = rep(0, 2))
    return(c(apply(sum_over_y, MARGIN = 1, sum), (1+x/beta)^{-alpha}))},
    FUN.VALUE = rep(0, 3))

  final_sum <- apply(sum_over_x, MARGIN = 1, sum)
  final_sum <- final_sum * delta
  final_sum[1:2] <- final_sum[1:2] * delta

  res <- final_sum[1]
  res_0 <- final_sum[2]
  first_mom_sq <- final_sum[3]^2

  if(cov){
    return(res-first_mom_sq)
  }else{
    return((res-first_mom_sq)/(res_0-first_mom_sq))
  }
}

#' @export
acf_trawl_revised <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp', cov=F){
  # Compute ACF with trawl process as latent
  seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
  trawl_fct <- GetTrawlFunctions(type)
  B1_func <- trawl_fct[[1]]
  B2_func <- trawl_fct[[2]]
  B3_func <- trawl_fct[[3]]

  b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))

  res_0 <- SquareMoment(xs = seq_kappa, delta = delta, beta = beta, b_oh = b_0_h, b_o_exc_h = b_0_minus_h)
  res <- CrossMoment(xs = seq_kappa, delta = delta, beta = beta, b_oh = b_0_h, b_o_exc_h = b_0_minus_h)
  first_mom_sq <- FirstMoment(xs = seq_kappa, delta = delta, beta = beta, b_oh = b_0_h, b_o_exc_h = b_0_minus_h)^2

  if(cov){
    return(res-first_mom_sq)
  }else{
    return((res-first_mom_sq)/(res_0-first_mom_sq))
  }
}

crossmoment_trawls <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp'){
  seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
  trawl_fct <- GetTrawlFunctions(type)
  B1_func <- trawl_fct[[1]]
  B2_func <- trawl_fct[[2]]
  B3_func <- trawl_fct[[3]]

  b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))

  res <- 0
  first_mom <- 0
  res_0 <- 0
  beta <- beta

  sum_over_x <- vapply(seq_kappa, FUN = function(x){
    x <- x + delta / 2
    sum_over_y <- vapply(X = seq_kappa, FUN = function(y){
      y <- y + delta / 2
      tmp1 <- (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
      return(tmp1)
    },
    FUN.VALUE = rep(0, 1))
    return(sum(sum_over_y))},
    FUN.VALUE = rep(0, 1))

  final_sum <- sum(sum_over_x)
  res <- final_sum * delta^2
  return(res)
}

acf_trawl_num_approx <- function(h, alpha, beta, kappa, rho, delta=0.5, type='exp', cov=T){
  vapply(h, function(h){
    acf_trawl(h, alpha = alpha, beta = beta, kappa = kappa,
              rho = rho, delta = delta, type = type, cov = cov)}, 1)
}

# acf_trawl_inv <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp', cov=T){
#   seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
#   c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type)
#   b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
#   b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
#   b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
#   res <- 0
#   first_mom <- 0
#   res_0 <- 0
#   beta <- beta
#
#   sum_over_x <- vapply(seq_kappa, FUN = function(x){
#     x <- x + delta / 2
#     # first_mom <- first_mom + (1+x/beta)^{-alpha}*(1+y/beta)^{-alpha}
#     sum_over_y <- vapply(X = seq_kappa, FUN = function(y){
#       y <- y + delta / 2
#       tmp_11 <- (1+(x-kappa)/beta)^{b_h_minus_0} * (1+(x+y-2*kappa)/beta)^{b_h_minus_0} * (1+(y-kappa)/beta)^{b_h_minus_0}
#       tmp_1e <- 2*(1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+(y-kappa)/beta)^{b_0_h}
#       tmp_ee <- (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
#       return(c(tmp_11, tmp_1e, tmp_ee))
#     },
#     FUN.VALUE = rep(0, 3))
#     return(c(apply(sum_over_y, MARGIN = 1, sum), (1+(x)/beta)^{-alpha}))},
#     FUN.VALUE = rep(0, 4))
#
#   # prob_xt_and_xs_pos <- 1 + (1+kappa/beta)^{b_0_minus_h}*(1+2*kappa/beta)^{b_0_h}*(1+kappa/beta)^{b_h_minus_0}
#   # prob_xt_and_xs_pos <- prob_xt_and_xs_pos - 2 * (1+kappa/beta)^{-alpha}
#
#   final_sum <- apply(sum_over_x, MARGIN = 1, sum)
#
#   res <- (final_sum[1] - final_sum[2] + final_sum[3]) * delta^2
#   res_0 <- sum(sum_over_x[4,]^2) * delta
#
#   # first_mom_sq <-  ((1+kappa/beta)^{-alpha}*(beta+kappa)/(alpha - 1))^2
#   first_mom_sq <- final_sum[4]^2 * delta^2
#   print(res)
#   print(res_0)
#   print(first_mom_sq)
#   if(cov){
#     return(res-first_mom_sq)
#   }else{
#     return((res-first_mom_sq)/(res_0-first_mom_sq))
#   }
# }
#
# acf_trawl_num_approx_inv <- function(h, alpha, beta, kappa, rho, delta=0.5, type='exp', cov=T){
#   vapply(h, function(h){
#     acf_trawl_inv(h, alpha = alpha, beta = beta, kappa = kappa,
#                   rho = rho, delta = delta, type = type, cov = cov)}, 1)}
