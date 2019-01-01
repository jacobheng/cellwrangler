#'Plot pattern correlation dendrogram
#'
#' plot_pattern_cor_dendro() plots a dendrogram showing the clustering of the correlation of a CoGAPS pattern 
#' set with a grouping variable
#'
#' @param CoGAPS_res_set a CoGAPS result set. If set to NULL, a Pmeans matrix may be supplied directly as input
#' for the Pattern_set parameter.
#' @param Pattern_set a string specifying the name of a pattern set (e.g. "nP20") in a CoGAPS_res_set. 
#' Alternatively, if CoGAPS_res_set is set to NULL, a Pmeans matrix may be supplied directly.
#' @param pattern_subset a vector specifying a subset of the pattern set to be plotted. Defaults to NULL.
#' @param annotation an annotation dataframe. Must contain columns with column names corresponding to values
#' in the group_vector. 
#' @param group_vector a vector of the grouping variable with values corresponding to column names in the
#' annotation object. Only required if cluster is set to "groups".
#' @param order_subgroups logical; if subgroups should be ordered according to clustering of first subgroup. 
#' Subgroups can be specified in a string with separation by " - ". Only valid if cluster is set to "groups".
#' @param cluster what to cluster. Possible values are "patterns" or "groups".
#' @param clustering_distance distance measure used in clustering. Possible values are "correlation",
#' "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski" (defined in dist() function.)
#' @param clustering_method clustering method used. Accepts the same values as hclust: "ward.D", 
#' "ward.D2", "single", "complete", "average" (=UPGMA), "mcquitty" (=WPGMA), "median" (=WPGMC) or 
#' "centroid" (=UPGMC).
#' 
#' @keywords plot_pattern_cor_dendro
#' @export
#' @return a ggplot object
#' @examples
#' plot_pattern_cor_heatmap(CoGAPS_res_set= my_CoGAPS_res, Pattern_set="nP30",  annotation=pData(cds), 
#' group_vector=group_vector, order_subgroups=T, cluster= "groups",
#' clustering_distance = "correlation", clustering_method = "complete")


#plot_pattern_cor_dendro
plot_pattern_cor_dendrogram <-  function(CoGAPS_res_set=NULL, Pattern_set, pattern_subset= NULL,   
                                         annotation, group_vector, order_subgroups=F, cluster="patterns", 
                                         clustering_distance="correlation", clustering_method="complete") {
  
  if(is.null(CoGAPS_res_set)==F) {
    Pmeans <- (t(CoGAPS_res_set[[Pattern_set]]$Pmean))} else { Pmeans <- Pattern_set }
  if(is.null(pattern_subset) == F) {
    Pmeans <- Pmeans[,pattern_subset]
  } else { Pmeans <- Pmeans }
  Pmeans <- Pmeans[unique(rownames(Pmeans)),]
  Pmeans_annot <- merge_by_rownames(annotation, Pmeans, all.x = F, all.y=T)
  for(i in c( (ncol(annotation)+1):(ncol(Pmeans)+ncol(annotation)) )){
    Pmeans_annot[,i]<-as.numeric(Pmeans_annot[,i])
  }
  patterns <- colnames(Pmeans_annot)[grep("Patt",colnames(Pmeans_annot))]
  
  #Calculate correlation
  tmp <- cor(Pmeans_annot[,patterns], Pmeans_annot[,group_vector])
  
  #cluster
  if(cluster == "patterns") {
    patterns.dendro <- as.dendrogram(pheatmap:::cluster_mat(tmp, distance=clustering_distance,method=clustering_method))
    d <- ggdendrogram(patterns.dendro, rotate=F)
    d.order <- order.dendrogram(patterns.dendro)
    
  } else { if(cluster=="groups") {
    ref_table <- as.data.frame(group_vector)
    colnames(ref_table) <- ("group_vector")
    #order_subgroups
    if(order_subgroups == T) {
      n_subgroups <- length(stringr::str_split(group_vector, pattern = " - ")[[1]])
      for(i in 1:n_subgroups) {
        ref_table[paste("subgroup",i,sep="")] <- as.factor(stringr::str_split_fixed(group_vector, pattern=" - ",
                                                                                    n=n_subgroups)[,i])
      }
      
      groups.dendro <- as.dendrogram(pheatmap:::cluster_mat(t(tmp[,unique(ref_table$subgroup1)]), 
                                                            distance=clustering_distance,method=clustering_method))
    
    } else {
      groups.dendro <- as.dendrogram(pheatmap:::cluster_mat(t(tmp), distance=cluster_groups_distance,
                                                            method=cluster_groups_method))
    }
    
    d <- ggdendrogram(groups.dendro, rotate=T) 
    d.order <- order.dendrogram(groups.dendro)
    
  }   else { print("Specify patterns or celltypes to cluster!") }
  }
  d.object <- list(d, d.order)
  return(d.object)
  
}