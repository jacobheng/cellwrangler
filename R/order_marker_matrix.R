#' Order expression matrix of marker genes by cell group
#'
#' @description order_marker_matrix() returns an expression matrix of marker genes ordered by cell group
#' 
#' @param cds CellDataSet object
#' @param marker_genes dataframe containing 2 columns: (1) gene symbol (2) cell group e.g. CellType
#' @param group_vector a vector classifying cells e.g. genotype or condition
#' @param scale "rows" or "cols" or "none"; indicates whether values should be centered and scaled in the 
#' row direction or the column direction or none.
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
#' @keywords order_marker_matrix()
#' @export
#' @return A matrix
#' @examples
#' order_marker_matrix(dat, marker_genes = myMarkers, group_vector = "CellType")


order_marker_matrix <- function(cds, marker_genes, group_vector, scale=T, cluster_rows = T, 
                                        cluster_cols = T, clustering_distance_rows = "euclidean", 
                                        clustering_method_rows = "complete", clustering_distance_cols = "euclidean", 
                                        clustering_method_cols = "complete") {
  
  marker_genes <- as.data.frame(marker_genes)  
  colnames(marker_genes) <- c("gene_symbol", "group")
  marker_genes$gene_id <- findGeneID(marker_genes$gene_symbol, cds)
  
  pData(cds)$group <- as.factor(pData(cds)[,group_vector])
  group_order <- as.data.frame(cbind(levels(pData(cds)$group), seq(1,length(levels(pData(cds)$group)), 1)))
  colnames(group_order) <- c("group","group_order")
  pData(cds)$group_order <- group_order$group_order[match(pData(cds)$group, group_order$group)]
  marker_genes$group_order <- group_order$group_order[match(marker_genes$group, group_order$group)]
  marker_genes$group_order <- as.numeric(as.character(marker_genes$group_order))
  
  
  heatData <-log10(exprs(cds)[marker_genes$gene_id,]+1)
  colnames(heatData)<-colnames(exprs(cds)[marker_genes$gene_id,])
  rownames(heatData)<-rownames(exprs(cds)[marker_genes$gene_id,])
  
  if(scale == T) {
    tmp <-as.matrix(t(apply(heatData,1,scale)))
    colnames(tmp)<-colnames(exprs(cds)[marker_genes$gene_id,])
    rownames(tmp)<-rownames(exprs(cds)[marker_genes$gene_id,])
  } else { tmp <- heatData }
  
  
  if (cluster_rows == T) {
    row.order <- order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(tmp, distance = clustering_distance_rows, 
                                                                       method = clustering_method_rows)))
    marker_genes$row_order <- as.numeric(as.character(row.order))
    tmp<-tmp[marker_genes[with(marker_genes, order(group_order,row_order)),]$gene_id,]
  }
  else {
    tmp<-tmp[marker_genes[with(marker_genes, order(group_order,row_order)),]$gene_id,]
  }
  
  if(cluster_cols == T) {
    col.order <- order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(t(tmp), distance = clustering_distance_cols, 
                                                                       method = clustering_method_cols)))
    pData(cds)$col_order <- col.order
    tmp<-tmp[,with(pData(cds),order(group_order,col_order))]
  }
  else {
    tmp<-tmp[,with(pData(cds),order(group_order))]
  }
  
  tmp<-tmp[marker_genes[with(marker_genes, order(group_order,row_order)),]$gene_id,]
  
  
  return(tmp)
  
}