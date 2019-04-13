#'Plot a clustered heatmap
#'
#' @description Draws a clustered heatmap using clustering methods implemented in the popular R package pheatmap. 
#' Returns a ggplot2 object
#' @param matrix a matrix of values
#' @param scale "rows" or "cols" or "none"; indicates whether values should be centered and scaled in the 
#' row direction or the column direction or none.
#' @param limits minimum and maximum values to define the range of the color scale
#' @param color_scale_type "diverge", "seq" or "rescale".
#' @param color_scale vector of colors to use for color scale; different options for color_scale_type requires
#' different number of colors - "diverge" (3 colors), "seq" (2 colors), "rescale" (3 colors)
#' @param scale_midpoint numeric; value of midpoint for "diverge" and "rescale" color scales
#' @param cluster_rows logical; if rows should be clustered
#' @param cluster_cols logical; if columns should be clustered
#' @param clustering_distance_rows distance measure used in clustering rows. Possible values are "correlation",
#' "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski" (defined in dist() function.)
#' @param clustering_method_rows clustering method used for rows. Accepts the same values as hclust: "ward.D", 
#' "ward.D2", "single", "complete", "average" (=UPGMA), "mcquitty" (=WPGMA), "median" (=WPGMC) or 
#' "centroid" (=UPGMC).
#' @param clustering_distance_cols distance measure used in clustering columns. Possible values are the same as
#' clustering_distance_rows.
#' @param clustering_method_cols clustering method used for columns. Accepts the same values as 
#' clustering_method_rows.
#' @param square_tiles logical; if tiles should be made square
#' @keywords plot_matrix_heatmap
#' @export
#' @return a ggplot2 object
#' @examples
#' plot_matrix_heatmap(matrix)


plot_matrix_heatmap <- function(matrix, scale="row", limits=NULL, color_scale_type = "diverge", 
                                color_scale=c("blue","lightyellow","red"), scale_midpoint=0, cluster_rows = T, 
                                cluster_cols =T, clustering_distance_rows="euclidean", 
                                clustering_method_rows="complete", clustering_distance_cols="euclidean", 
                                clustering_method_cols="complete", square_tiles=F) {
  tmp <- matrix
  if(scale=="rows"){
    tmp <- t(apply(tmp,1,scale))
    colnames(tmp) <- colnames(matrix)
  } else { if(scale=="cols"){  
    tmp <- (apply(tmp,2,scale))
    colnames(tmp) <- colnames(matrix)
  } else{ tmp <- as.matrix(tmp) } } 
  matrix_melt <- reshape2:::melt.matrix(tmp)
  colnames(matrix_melt) <- c("row.name","sample","value")
  if(cluster_rows == T) {
    row.order<-order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(tmp,distance=clustering_distance_rows,method=clustering_method_rows)))
    matrix_melt$row.name <-factor(matrix_melt$row.name,levels=rev(matrix_melt$row.name[row.order]))
  } else { matrix_melt$row.name<-factor(matrix_melt$row.name,levels=rownames(matrix))  }
  if(cluster_cols == T) {
    col.order<-order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(t(tmp), distance=clustering_distance_cols,method=clustering_method_cols)))
    matrix_melt$sample<-factor(matrix_melt$sample,levels=colnames(matrix)[col.order])
  } else { matrix_melt$sample<-factor(matrix_melt$sample,levels=colnames(matrix)) }
  p<-ggplot(matrix_melt,aes(x=sample,y=row.name,fill=value))
  p <- p + geom_tile(size=0.25) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
  if(square_tiles == TRUE) {
    p <- p + coord_equal()
  } else { p <- p }
  if(color_scale_type=="diverge") {
    p <- p + scale_fill_gradient2(low=color_scale[1],mid=color_scale[2],high=color_scale[3], midpoint= scale_midpoint, limits=limits, oob=scales::squish)
  } else {if(color_scale_type=="seq") {
    scale_midpoint
    p <- p + scale_fill_gradient(low=color_scale[1],high=color_scale[2],limits=limits,oob=scales::squish)
  }else {if(color_scale_type=="rescale") {
    print("rescale")
    p <- p + scale_fill_gradientn(colours =  color_scale, values=scales::rescale(c(min(matrix_melt$value),scale_midpoint, max(matrix_melt$value))))
  } else { p <- p } } }
  return(p)
}