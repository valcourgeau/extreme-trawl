autocovariance_matrix <- function(score, params, k) {
  assertthat::assert_that(dim(score)[1], 4)
  return(acf(score, lag.max = k, type = "covariance")$acf)
}

make_hac <- function(acf_matrices, near_pd = F) {
  # acf_matrices is k x n_params x n_params
  k <- dim(acf_matrices)[1]
  weights <- 1 - (seq_len(k) - 1) / (k + 1)
  weighted_acf_matrices <- lapply(
    seq_len(k), function(i) {
      weights[i] * acf_matrices[i, , ]
    }
  )
  hac_ests <- Reduce(`+`, weighted_acf_matrices)
  if (near_pd) {
    hac_ests <- as.matrix(Matrix::nearPD(hac_ests))
  }
  return(hac_ests)
}
