autocovariance_matrix <- function(score, k, params_indices = NA) {
  if (is.na(params_indices)) {
    if (is.matrix(score)) {
      params_indices <- seq_len(dim(score)[2])
    } else {
      params_indices <- seq_len(1)
    }
  }
  acf_vals <- acf(score, lag.max = k, type = "covariance")$acf
  acf_vals <- acf_vals[, params_indices, params_indices]
  return(acf_vals)
}

make_hac <- function(acf_matrices, near_pd = F) {
  # acf_matrices is k x n_params x n_params
  k <- dim(acf_matrices)[1]
  weights <- 1 - (seq_len(k) - 1) / (k + 1)
  weighted_acf_matrices <- lapply(
    seq_len(k), function(i) {
      tmp <- acf_matrices[i, , ]
      if (i > 1) {
        tmp <- weights[i] * (tmp + t(tmp))
      }
      return(tmp)
    }
  )
  hac_ests <- Reduce(`+`, weighted_acf_matrices)
  assertthat::assert_that(all(dim(hac_ests) == dim(acf_matrices)[c(2, 3)]))
  if (near_pd) {
    hac_ests <- as.matrix(Matrix::nearPD(hac_ests))
  }
  return(hac_ests)
}
