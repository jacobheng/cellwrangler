#'Sparse prcomp irlba implemented in Monocle 3 alpha
#'
#' @description Efficient computation of a truncated principal components analysis of a given data matrix
#' using an implicitly restarted Lanczos method from the \code{\link{irlba}} package.The augmented implicitly restarted Lanczos bidiagonalization algorithm (IRLBA) finds a few 
#' approximate largest (or, optionally, smallest) singular values and corresponding singular vectors of a 
#' sparse or dense matrix using a method of Baglama and Reichel. It is a fast and memory-efficient way to 
#' compute a partial SVD.
#' @param x a numeric or complex matrix (or data frame) which provides
#'          the data for the principal components analysis.
#' @param retx a logical value indicating whether the rotated variables should be returned.
#' @param center a logical value indicating whether the variables should be
#'          shifted to be zero centered. Alternately, a centering vector of length
#'          equal the number of columns of \code{x} can be supplied.
#' @param scale. a logical value indicating whether the variables should be
#'          scaled to have unit variance before the analysis takes place.
#'          The default is \code{FALSE} for consistency with S, but scaling is often advisable.
#'          Alternatively, a vector of length equal the number of columns of \code{x} can be supplied.
#'
#'          The value of \code{scale} determines how column scaling is performed
#'          (after centering).  If \code{scale} is a numeric vector with length
#'          equal to the number of columns of \code{x}, then each column of \code{x} is
#'          divided by the corresponding value from \code{scale}.  If \code{scale} is
#'          \code{TRUE} then scaling is done by dividing the (centered) columns of
#'          \code{x} by their standard deviations if \code{center=TRUE}, and the
#'          root mean square otherwise.  If \code{scale} is \code{FALSE}, no scaling is done.
#'          See \code{\link{scale}} for more details.
#' @param n integer number of principal component vectors to return, must be less than
#' \code{min(dim(x))}.
#' @param ... additional arguments passed to \code{\link{irlba}}.

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