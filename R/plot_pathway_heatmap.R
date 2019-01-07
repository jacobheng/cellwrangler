#'Plot pathway heatmap
#'
#'plot_pathway_heatmap() plots the results of output from the fgsea() function.
#' 
#' @param path_table a dataframe of values
#' @param path_res_df output from the fgsea() function
#' @param color_scale vector of colors to use for color scale; different options for color_scale_type requires
#' different number of colors - "diverge" (3 colors), "seq" (2 colors), "rescale" (3 colors)
#' @param cluster_rows logical; if rows should be clustered
#' @param cluster_cols logical; if columns should be clustered
#' @param mark_sig logical; whether to mark statistically significant results with a "*".
#' @param sig_para parameter for statistical significance.
#' @param sigCutoff cutoff for statistical significance.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible values are "correlation",
#' "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski" (defined in dist() function.)
#' @param clustering_method_rows clustering method used for rows. Accepts the same values as hclust: "ward.D", 
#' "ward.D2", "single", "complete", "average" (=UPGMA), "mcquitty" (=WPGMA), "median" (=WPGMC) or 
#' "centroid" (=UPGMC).
#' @param clustering_distance_cols distance measure used in clustering columns. Possible values are the same as
#' clustering_distance_rows.
#' @param clustering_method_cols clustering method used for columns. Accepts the same values as 
#' clustering_method_rows.
#' @keywords plot_pathway_heatmap()
#' @export
#' @return A ggplot2 object
#' @examples
#' plot_pathway_heatmap(path_table=NULL, path_res_df=my_pathway_enrichment)


plot_pathway_heatmap <- function(path_table=NULL, path_res_df, color_scale=c("blue","lightyellow","red"), 
                                 cluster_row = T, cluster_col =T , mark_sig = T, sig_para="pval",sigCutoff=0.05, 
                                 clustering_distance_rows="euclidean",clustering_method_rows="complete", 
                                 clustering_distance_cols="euclidean", clustering_method_cols="complete") {
  if(is.null(path_table)==TRUE){
    path_table <- dcast(path_res_df[,c("pathway","NES","CellType")], pathway~CellType, value.var = "NES")
    rownames(path_table) <- path_table$pathway
    path_table <- path_table[,-1]
    path_table[is.na(path_table)] <- 0 } else {
      path_table <- path_table
    }
  if(cluster_row == T) {
    row.order<-order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(path_table,distance=clustering_distance_rows,method=clustering_method_rows)))
    path_res_df$pathway<-factor(path_res_df$pathway,levels=rownames(path_table)[row.order])
  } else { path_res_df <- path_res_df }
  if(cluster_col == T) {
    col.order<-order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(t(path_table), distance=clustering_distance_cols,method=clustering_method_cols)))
    path_res_df$CellType<-factor(path_res_df$CellType,levels=colnames(path_table)[col.order])
  } else { path_res_df <- path_res_df }
  path_res_df$sig <- ifelse(path_res_df[sig_para]<=sigCutoff,"*","")
  p<-ggplot(path_res_df,aes(x=CellType,y=pathway,fill=NES))
  p <- p + geom_tile(size=0.25) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + scale_fill_gradient2(low=color_scale[1],mid=color_scale[2],high=color_scale[3],na.value = "grey50")
  if(mark_sig == TRUE) {
    p<-p + geom_text(aes(label=sig),nudge_x=0)
  } else {
    p<-p
  }
  return(p)
}