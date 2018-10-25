



sparse_pca <- function (x, n_pcs, mu = NULL, s = NULL, center_scale = TRUE) 
{
  if (is.null(mu) && center_scale) 
    mu <- colMeans(x)
  if (is.null(s) && center_scale) 
    s <- apply(x, 2, sd, na.rm = TRUE)
  irlba_wrapper <- function(...) {
    res <- try(irlba::irlba(fastpath = FALSE, ...), silent = T)
    if (inherits(res, "try-error")) {
      res <- irlba::irlba(...)
    }
    res
  }
  if (center_scale) {
    s[s == 0] <- min(s[s > 0])
    svd_res <- irlba_wrapper(x, n_pcs, center = mu, scale = s)
  }
  else {
    svd_res <- irlba_wrapper(x, n_pcs)
  }
  n <- dim(x)[1]
  variance_sum <- sum(apply(x, 2, var, na.rm = TRUE)/(s^2))
  var_pcs <- svd_res$d^2/(n - 1)/variance_sum
  return(list(x = svd_res$u %*% diag(svd_res$d), rotation = svd_res$v, 
              sdev = svd_res$d/sqrt(n - 1), tot_var = variance_sum, 
              var_pcs = var_pcs))
}