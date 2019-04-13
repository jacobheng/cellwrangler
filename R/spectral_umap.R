
#'Reduce dimensions of a matrix with pca and UMAP
#'
#' @description spectral_umap() performs dimensionality reduction using the package monocle wrapper function around
#' the python implementation of UMAP (Uniform Manifold Approximation and Projection). If no prcomp_object 
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
#' @param umap_version UMAP implementations to use; options are "default", "monocle" or "uwot". "monocle" only
#' works if monocle 3 alpha and above is installed. The default option uses the UMAP function implemented
#' in monocle 3 alpha, and works even without monocle 3 alpha installed.
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

spectral_umap <- function(matrix, log_matrix=TRUE, prcomp_object=NULL, pca_version="default", center=T,scale=T, 
                          dims=1:10, umap_version="default", n_neighbors = 30L, metric= "correlation", min_dist = 0.1, 
                          spread = 1) {
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
        
        pca_res <- cellrwangler::monocle_sparse_prcomp_irlba(Matrix::t(tmp), n = min(max(dims), min(dim(tmp)) - 1), 
                                       center = center, scale. = scale)
        
      } else { print("need to specify pca version!")}
    }
  }
  
  
  if(umap_version=="default") {
  umap_proj <- cellwrangler::monocle_UMAP(X=pca_res$x[,dims], log=F, n_neighbors = n_neighbors, metric = metric,
                             min_dist = min_dist, spread = spread)
  colnames(umap_proj) <- c("UMAP.1", "UMAP.2")
  #rownames(umap_proj) <- rownames(matrix)
  } else {
    if(umap_version=="monocle") {
      umap_proj <- monocle::UMAP(X=pca_res$x[,dims], log=F, n_neighbors = n_neighbors, metric = metric,
                        min_dist = min_dist, spread = spread)
      colnames(umap_proj) <- c("UMAP.1", "UMAP.2")
      #rownames(umap_proj) <- rownames(matrix)
    } else {
      if(umap_version=="uwot") {
        umap_proj <- uwot::umap(pca_res$x[,dims], n_neighbors = n_neighbors, metric = metric,
                                   min_dist = min_dist, spread = spread)
        colnames(umap_proj) <- c("UMAP.1", "UMAP.2")
        #rownames(umap_proj) <- rownames(matrix)
      } else { print("Must specify implementation as default, monocle or uwot!")}
    }
  }
  return(umap_proj)
} 

