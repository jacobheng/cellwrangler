#'Plot a clustered dendrogram
#'
#' @description Draws a clustered heatmap using clustering methods implemented in the popular R package pheatmap. 
#' Returns a ggplot2 object
#' @param matrix a matrix of values
#' @param scale "rows" or "cols" or "none"; indicates whether values should be centered and scaled in the 
#' row direction or the column direction or none.
#' @param cluster "rows" or "cols"; whether to cluster rows or columns.
#' @param clustering_distance distance measure to be used. Possible values are "correlation",
#' "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski" (defined in dist() function.)
#' @param clustering_method clustering method to be used. Accepts the same values as hclust: "ward.D", 
#' "ward.D2", "single", "complete", "average" (=UPGMA), "mcquitty" (=WPGMA), "median" (=WPGMC) or 
#' "centroid" (=UPGMC).
#' @keywords plot_matrix_heatmap
#' @export
#' @return a ggplot2 object
#' @examples
#' plot_matrix_heatmap(matrix)

plot_matrix_dendrogram <-  function(matrix, scale="row",cluster="col", clustering_distance="euclidean", 
                                    clustering_method="complete") {
  tmp <- matrix
  if(scale=="rows"){
    tmp <- t(apply(tmp,1,scale))
    colnames(tmp) <- colnames(matrix)
  } else { if(scale=="cols"){  
    tmp <- (apply(tmp,2,scale))
    colnames(tmp) <- colnames(matrix)
  } else{ tmp <- as.matrix(tmp) } } 
  if(cluster == "cols") {
    col.dendro <- stats::as.dendrogram(pheatmap:::cluster_mat(t(tmp), distance=clustering_distance,
                                                              method=clustering_method))
    d <- ggdendrogram(col.dendro, rotate=F) 
  } else { if(cluster=="rows") {
    row.dendro <- stats::as.dendrogram(pheatmap:::cluster_mat(tmp, distance=clustering_distance,
                                                              method=clustering_method))
    d <- ggdendrogram(row.dendro, rotate=T) 
  }  else { print("Specify col or row to cluster!") }
  }
  return(d)
}