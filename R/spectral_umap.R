
#'Reduce dimensions of a matrix with pca and UMAP
#'
#' spectral_umap() performs dimensionality reduction using the package monocle wrapper function around
#' the python implementation of UMAP (Uniform Manifold Approximation and Projection). If no prcomp_object 
#' is supplied, principal component analysis will be performed using the irlba package.
#' @param matrix a matrix of values to perform dimensionality reduction on; by default, rows are genes
#' and columns are cells
#' @param log_matrix if log10 transformation is to be performed on the matrix; defaults to TRUE
#' @param prcomp_obj a principal component analysis object produced by the prcomp or prcomp_irlba functions;
#' if no object is supplied, irlba_prcomp will be run on the matrix to return 50 dimensions; defaults to NULL; 
#' if a prcomp object is supplied, matrix is not required
#' @param dims dimensions from the prinicpal component analysis to use; defaults to 1:10 (i.e. 1st to 10th
#' principal components)
#' @param n_neighbors float (optional, default 15) The size of local neighborhood (in terms of number of 
#' neighboring sample points) used for manifold approximation. Larger values result in more global views of 
#' the manifold, while smaller values result in more local data being preserved. In general values should be 
#' in the range 2 to 100.
#' @param metric string or function (optional, default 'correlation') The metric to use to compute distances in 
#' high dimensional space. If a string is passed it must match a valid predefined metric. If a general metric is 
#' required a function that takes two 1d arrays and returns a float can be provided. For performance purposes it 
#' is required that this be a numba jit'd function. Valid string metrics include: * euclidean * manhattan * 
#' chebyshev * minkowski * canberra * braycurtis * mahalanobis * wminkowski * seuclidean * cosine * correlation 
#' * haversine * hamming * jaccard * dice * russelrao * kulsinski * rogerstanimoto * sokalmichener * sokalsneath 
#' * yule Metrics that take arguments (such as minkowski, mahalanobis etc.) can have arguments passed via the 
#' metric_kwds dictionary. At this time care must be taken and dictionary elements must be ordered appropriately; 
#' this will hopefully be fixed in the future.
#' @param min_dist float (optional, default 0.1) The effective minimum distance between embedded points. 
#' Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are 
#' drawn closer together, while larger values will result on a more even dispersal of points. The value should 
#' be set relative to the “spread“ value, which determines the scale at which embedded points will be spread out.
#' @param spread float (optional, default 1.0) The effective scale of embedded points. In combination with 
#' “min_dist“ this determines how clustered/clumped the embedded points are.
#' @keywords spectral_umap
#' @export
#' @return A matrix with two columns containing coordinates of each row for two dimensions respectively
#' @examples
#' spectral_umap(matrix, log_matrix=TRUE, prcomp_object=NULL, dims=1:10)

spectral_umap <- function(matrix, log_matrix=TRUE, prcomp_object=NULL, dims=1:10, n_neighbors = 30L, 
                          metric= "correlation", min_dist = 0.1, spread = 1) {
  if(is.null(prcomp_object) == FALSE) {
    pca_res <- prcomp_object
  } else {
    tmp <- matrix
    if(log_matrix==TRUE) {
      tmp <- log10(matrix+1)
    } else { tmp <- tmp}
    pca_res <- cellrangerRkit:::sparse_pca(t(tmp), n_pcs=40, center_scale = T) 
    }
  umap_proj <- monocle::UMAP(pca_res$x[,dims], log=F, n_neighbors = n_neighbors, metric = metric,
                             min_dist = min_dist, spread = spread)
  colnames(umap_proj) <- c("UMAP.1", "UMAP.2")
  return(umap_proj)
} 