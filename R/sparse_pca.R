#'Maintain sparsity and run a fast PCA algorithm
#'
#' sparse_pca() runs a fast PCA algorithm on a sparse matrix. Adapted from CellRangerRkit.
#'
#' @param matrix a sparse matrix
#' @param n_pcs number of principal components to compute
#' @param mu column means
#' @param s column standard deviations
#' @param center_scale logical; whether to perform centering and scaling
#' 
#' @keywords sparse_pca()
#' @export
#' @return an irlba pca object
#' @examples
#' pca_object <- sparse_pca(exprs(dat), n_pcs = 40)



sparse_pca <- function (matrix, n_pcs, mu = NULL, s = NULL, center_scale = TRUE) 
{
  if (is.null(mu) && center_scale) 
    mu <- colMeans(matrix)
  if (is.null(s) && center_scale) 
    s <- apply(matrix, 2, sd, na.rm = TRUE)
  irlba_wrapper <- function(...) {
    res <- try(irlba::irlba(fastpath = FALSE, ...), silent = T)
    if (inherits(res, "try-error")) {
      res <- irlba::irlba(...)
    }
    res
  }
  if (center_scale) {
    s[s == 0] <- min(s[s > 0])
    svd_res <- irlba_wrapper(matrix, n_pcs, center = mu, scale = s)
  }
  else {
    svd_res <- irlba_wrapper(matrix, n_pcs)
  }
  n <- dim(matrix)[1]
  variance_sum <- sum(apply(matrix, 2, var, na.rm = TRUE)/(s^2))
  var_pcs <- svd_res$d^2/(n - 1)/variance_sum
  return(list(matrix = svd_res$u %*% diag(svd_res$d), rotation = svd_res$v, 
              sdev = svd_res$d/sqrt(n - 1), tot_var = variance_sum, 
              var_pcs = var_pcs))
}