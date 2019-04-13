#'Reduce dimensions of a matrix with pca and tSNE
#'
#' @description spectral_tsne() performs dimensionality reduction using the package Rtsne's wrapper for the 
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
        
        pca_res <- cellwrangler::monocle_sparse_prcomp_irlba(Matrix::t(tmp), n = min(max(dims), min(dim(tmp)) - 1), 
                                               center = center, scale. = scale)
        
      } else { print("need to specify pca version!")}
    }
    }
  Rtsne_obj <- Rtsne::Rtsne(pca_res$x[,dims],perplexity=perplexity, pca = FALSE)
  tsne_proj <- Rtsne_obj$Y
  colnames(tsne_proj) <- c("TSNE.1", "TSNE.2")
  return(tsne_proj)
} 