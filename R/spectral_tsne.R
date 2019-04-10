#'Reduce dimensions of a matrix with pca and tSNE
#'
#' spectral_tsne() performs dimensionality reduction using the package Rtsne's wrapper for the 
#' C++ implementation of Barnes-Hut t-Distributed Stochastic Neighbor Embedding. If no prcomp_object 
#' is supplied, principal component analysis will be performed using an irlba algorithm for sparse data.
#' @param matrix a matrix of values to perform dimensionality reduction on; by default, rows are genes
#' and columns are cells
#' @param log_matrix if log10 transformation is to be performed on the matrix; defaults to TRUE
#' @param prcomp_obj a principal component analysis object produced by the prcomp or prcomp_irlba functions;
#' if no object is supplied, sparse_pca will be run on the matrix to return 50 dimensions; defaults to NULL;
#' if a prcomp object is supplied, matrix is not required
#' @param pca_version PCA implementation to use. Possible values are "default" for sparse_pca() or "monocle"
#' for the sparse_irlba_prcomp implemented in Monocle 3 alpha.
#' @param center a logical value indicating whether the variables should be shifted to be zero centered. 
#' Alternately, if the "monocle" pca version is used, a centering vector of length equal the number of columns of 
#' x can be supplied.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance 
#' before the analysis takes place. If the "default" pca version is used and center = TRUE, scaling will be
#' also default to TRUE. Alternatively, if the "monocle" pca version is used, a vector of length equal the 
#' number of columns of x can be supplied.
#' @param dims dimensions from the prinicpal component analysis to use; defaults to 1:10 (i.e. 1st to 10th
#' principal components)
#' @param perplexity numeric; perplexity parameter for tSNE; defaults to 30
#' @keywords spectral_tsne
#' @export
#' @return A matrix with two columns containing coordinates of each row for two dimensions respectively
#' @examples
#' spectral_tsne(matrix, log_matrix=TRUE, prcomp_object=NULL, dims=1:10, perplexity=30)

spectral_tsne <- function(matrix, log_matrix=TRUE, prcomp_object=NULL, dims=1:10, perplexity=30) {
  if(is.null(prcomp_object) == FALSE) {
    pca_res <- prcomp_object
   } else { 
    tmp <- matrix
    if(log_matrix==TRUE) {
      tmp <- log10(matrix+1)
    } else { tmp <- tmp}
    if(pca_version=="default") {
      pca_res <- cellwrangler:::sparse_pca(Matrix::t(tmp), n_pcs=max(dims), center_scale = center) 
    } else {
      if(pca_version=="monocle") {
        
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
        
        
        pca_res <- monocle_sparse_prcomp_irlba(Matrix::t(tmp), n = min(max(dims), min(dim(tmp)) - 1), 
                                               center = center, scale. = scale)
        
      } else { print("need to specify pca version!")}
    }
    }
  Rtsne_obj <- Rtsne::Rtsne(pca_res$x[,dims],perplexity=perplexity, pca = FALSE)
  tsne_proj <- Rtsne_obj$Y
  colnames(tsne_proj) <- c("TSNE.1", "TSNE.2")
  return(tsne_proj)
} 