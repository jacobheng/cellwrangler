#'Perform sparse pca on a matrix
#'
#' sparse_pca() performs a fast pca on a matrix while maintaing sparsity.
#' @param matrix a matrix of values to perform dimensionality reduction on; by default, rows are genes
#' and columns are cells
#' @param n_pcs number of prinicpal components to compute
#' @param mu column means
#' @param s column standard deviations
#' @param center_scale perform centering and scaling
#' @export
#' @return A list containing
#' "x" - The rotated data matrix where rows are barcodes and columns are PCs
#' "sdev" - the standard deviations of the principal components (i.e., sqrt of eigvals of the covariance matrix)
#' "rotation" - The loadings (eigenvectors) where each column is a PC
#' "tot_var" - The total variation in the scaled and centered matrix (this is also the effective rank of the matrix)
#' "var_pcs" - The proportion of variance explained by each principle comoponent
#' @examples
#' sparse_pca(mtx, n_pcs=10)

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