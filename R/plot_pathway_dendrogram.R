#'Plot pathway dendrogram
#'
#'@description plot_pathway_dendrogram() plots the results of output from the fgsea() function.
#' 
#' @param path_res_df dataframe should contain columns "CellType", "pathway", "NES";can be output from 
#' the fgsea() function
#' @param cluster what to cluster; possible values are "pathway" or "CellType"
#' @param clustering_distance distance measure used in clustering. Possible values are "correlation",
#' "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski" (defined in dist() function.)
#' @param clustering_method clustering method used for clustering. Accepts the same values as hclust: "ward.D", 
#' "ward.D2", "single", "complete", "average" (=UPGMA), "mcquitty" (=WPGMA), "median" (=WPGMC) or 
#' "centroid" (=UPGMC).
#' @param flip_dendro logical; if dendrogram should be flipped
#' @keywords plot_pathway_dendrogram()
#' @export
#' @return A ggplot2 object
#' @examples
#' plot_pathway_dendrogram(path_res_df=my_pathway_enrichment)

plot_pathway_dendrogram <-  function(path_res_df, cluster="col",clustering_distance="euclidean", 
                                     clustering_method="complete", flip_dendro = F) {
  path_table <- dcast(path_res_df[,c("pathway","NES","CellType")], pathway~CellType, value.var = "NES")
  rownames(path_table) <- path_table$pathway
  path_table <- path_table[,-1]
  path_table[is.na(path_table)] <- 0 
  if(cluster == "CellType") {
    col.dendro <- as.dendrogram(pheatmap:::cluster_mat(t(path_table), distance=clustering_distance, method=clustering_method))
    if(flip_dendro == T) {
    d <- ggdendrogram(col.dendro, rotate=T) 
    } else { d <- ggdendrogram(col.dendro, rotate=F) }
  } else { if(cluster=="pathway") {
    row.dendro <- as.dendrogram(pheatmap:::cluster_mat(path_table, distance=clustering_distance, method=clustering_method))
    if(flip_dendro == T) {
    d <- ggdendrogram(row.dendro, rotate=F) 
    } else { d <- ggdendrogram(row.dendro, rotate=T)  }
  }  else { print("Specify col or row to cluster!") }
  }
  return(d)
}