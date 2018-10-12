#'Reduce dimensions of a matrix with pca and tSNE
#'
#' spectral_tsne() performs dimensionality reduction using the package Rtsne's wrapper for the 
#' C++ implementation of Barnes-Hut t-Distributed Stochastic Neighbor Embedding. If no prcomp_object 
#' is supplied, principal component analysis will be performed using the irlba package.
#' @param matrix a matrix of values to perform dimensionality reduction on
#' @param log_matrix if log10 transformation is to be performed on the matrix; defaults to TRUE
#' @param prcomp_obj a principal component analysis object produced by the prcomp or prcomp_irlba functions;
#' if no object is supplied, irlba_prcomp will be run on the matrix to return 50 dimensions; defaults to NULL
#' @param dims dimensions from the prinicpal component analysis to use; defaults to 1:10 (i.e. 1st to 10th
#' principal components)
#' @param perplexity numeric; perplexity parameter for tSNE; defaults to 30
#' @keywords spectral_tsne
#' @export
#' @return A matrix with two columns containing coordinates of each row for two dimensions respectively
#' @examples
#' spectral_tsne(matrix, log_matrix=TRUE, prcomp_object=NULL, dims=1:10, perplexity=30)

spectral_tsne <- function(matrix, log_matrix=TRUE, prcomp_object=NULL, dims=1:10, perplexity=30) {
  tmp <- matrix
  if(log_matrix==TRUE) {
    tmp <- log10(matrix+1)
  } else { tmp <- tmp}
  if(is.null(prcomp_object) == FALSE) {
    pca_res <- prcomp_object
  } else { pca_res <- irlba::prcomp_irlba(t(tmp), n=50) }
  Rtsne_obj <- Rtsne::Rtsne(pca_res$x[,dims],perplexity=perplexity)
  tsne_proj <- Rtsne_obj$Y
  colnames(tsne_proj) <- c("TSNE.1", "TSNE.2")
  return(tsne_proj)
} 