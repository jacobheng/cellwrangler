




monocle_sparse_prcomp_irlba <- function (x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, 
          ...) 
{
  a <- names(as.list(match.call()))
  ans <- list(scale = scale.)
  if ("tol" %in% a) 
    warning("The `tol` truncation argument from `prcomp` is not supported by\n            `prcomp_irlba`. If specified, `tol` is passed to the `irlba` function to\n            control that algorithm's convergence tolerance. See `?prcomp_irlba` for help.")
  orig_x = x
  if (class(x) != "DelayedMatrix") 
    x = DelayedArray(x)
  args <- list(A = orig_x, nv = n)
  if (is.logical(center)) {
    if (center) 
      args$center <- DelayedMatrixStats::colMeans2(x)
  }
  else args$center <- center
  if (is.logical(scale.)) {
    if (is.numeric(args$center)) {
      scale. = sqrt(DelayedMatrixStats::colVars(x))
      if (ans$scale) 
        ans$totalvar <- ncol(x)
      else ans$totalvar <- sum(scale.^2)
    }
    else {
      if (ans$scale) {
        scale. = sqrt(DelayedMatrixStats::colSums2(x^2)/(max(1, 
                                                             nrow(x) - 1L)))
        ans$totalvar = sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 
                                                                                  1L)))
      }
      else {
        ans$totalvar = sum(DelayedMatrixStats::colSums2(x^2)/(nrow(x) - 
                                                                1L))
      }
    }
    if (ans$scale) 
      args$scale <- scale.
  }
  else {
    args$scale <- scale.
    ans$totalvar = sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 
                                                                              1L)))
  }
  if (!missing(...)) 
    args <- c(args, list(...))
  s <- do.call(irlba, args = args)
  ans$sdev <- s$d/sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), 
                                  sep = "")
  ans$center <- args$center
  if (retx) {
    ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN = `*`)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), 
                             sep = "")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}