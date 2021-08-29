
#' Autocovariance matrix from vector.
#' @param score Vector series.
#' @param k Maximum horizon.
#' @return Autocovariance matrix
#' @examples
#' d <- 5
#' n <- 1000
#' vals <- matrix(norm(n * d), ncol = d)
#' make_hac(vals, 5)
#' @export
#' @importFrom stats acf
autocovariance_matrix <- function(score, k) {
  return(acf(score, lag.max = k, type = "covariance")$acf)
}


#' Autocorrelation-proof estimates using HAC estimator.
#' @param acf_matrices Matrices to use in the estimators.
#' @param near_pd Boolean, nearest positive definite matrix from HAC estimate.
#' @return HAC-transformed matrix
#' @examples
#' d <- 5
#' mats <- array(0, c(d, d, d))
#' mats <- apply(mats, c(2, 3), function(x) {
#'   diag(x) <- 1
#'   return(x)
#' })
#' make_hac(mats)
#' @importFrom Matrix nearPD
#' @importFrom assertthat assert_that
#' @export
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
