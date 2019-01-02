
#'Plot pattern correlation heatmap
#'
#' plot_pattern_cor_heatmap() plots a heatmap showing the correlation of a CoGAPS pattern set with a grouping
#' variable
#'
#' @param CoGAPS_res_set a CoGAPS result set. If set to NULL, a Pmeans matrix may be supplied directly as input
#' for the Pattern_set parameter.
#' @param Pattern_set a string specifying the name of a pattern set (e.g. "nP20") in a CoGAPS_res_set. 
#' Alternatively, if CoGAPS_res_set is set to NULL, a Pmeans matrix may be supplied directly.
#' @param pattern_subset a vector specifying a subset of the pattern set to be plotted. Defaults to NULL.
#' @param annotation an annotation dataframe. Must contain a group_variable column with values corresponding to 
#' those in group_vector.
#' @param group_variable column name of column in annotation object with values corresponding to those 
#' in group_vector.
#' @param group_vector a vector of with pre-specified values corresponding to those in the group_variable in
#' the annotation object.
#' @param cluster_groups logical; if groups should be clustered.
#' @param order_subgroups logical; if subgroups should be ordered according to clustering of first subgroup. 
#' Subgroups can be specified in a string with separation by " - ".
#' @param cluster_groups_distance distance measure used in clustering groups. Possible values are "correlation",
#' "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski" (defined in dist() function.)
#' @param cluster_groups_method clustering method used for groups. Accepts the same values as hclust: "ward.D", 
#' "ward.D2", "single", "complete", "average" (=UPGMA), "mcquitty" (=WPGMA), "median" (=WPGMC) or 
#' "centroid" (=UPGMC).
#' @param cluster_patterns logical; if patterns should be clustered.
#' @param cluster_patterns_vector a vector specifying the order of patterns. If NULL, patterns will be clustered
#' according to distance and method specified below.
#' @param cluster_patterns_distance distance measure used in clustering patterns. Possible values are the
#' same as cluster_groups_distance().
#' @param cluster_patterns_method clustering method used for patterns. Possible values are the
#' same as cluster_groups_method().
#' @param limits limits of the color scale for the heatmap.
#' @param square_tiles logical; if tiles should be made square
#' 
#' @keywords plot_pattern_cor_heatmap
#' @export
#' @return a ggplot object
#' @examples
#' plot_pattern_cor_heatmap(CoGAPS_res_set = my_CoGAPS_res, Pattern_set = "nP30",annotation=pData(cds), 
#' group_variable="genotype", group_vector=genotype_vector, cluster_groups=T, order_subgroups=T, 
#' cluster_groups_distance = "correlation", cluster_groups_method = "complete", cluster_patterns= T, 
#' cluster_patterns_distance = "correlation", cluster_patterns_method = "complete")


#plot_pattern_cor_heatmap
plot_pattern_cor_heatmap <- function(CoGAPS_res_set=NULL, Pattern_set, pattern_subset= NULL, annotation, group_variable, 
                                     group_vector, cluster_groups=F, order_subgroups=F ,cluster_groups_distance="correlation", 
                                     cluster_groups_method="complete", cluster_patterns=F, 
                                     cluster_patterns_vector=NULL, cluster_patterns_distance="correlation", 
                                     cluster_patterns_method="complete",limits=NULL, square_tiles=T) {
  if(is.null(CoGAPS_res_set)==F) {
    Pmeans <- (t(CoGAPS_res_set[[Pattern_set]]@sampleFactors))} else { Pmeans <- Pattern_set }
  if(is.null(pattern_subset) == F) {
    Pmeans <- Pmeans[,pattern_subset]
  } else { Pmeans <- Pmeans }
  Pmeans <- Pmeans[unique(rownames(Pmeans)),]
  #Create columns in annotation according to group_vector
  lapply(group_vector, function(x){
    tmp <- annotation[group_variable] == x
    annotation[[x]] <<- as.numeric(tmp)
  })
  #Merge Pmeans and annotation
  Pmeans_annot <- merge_by_rownames(annotation, Pmeans, all.x = F, all.y=T)
  for(i in c( (ncol(annotation)+1):(ncol(Pmeans)+ncol(annotation)) )){
    Pmeans_annot[,i]<-as.numeric(Pmeans_annot[,i])
  }
  patterns <- colnames(Pmeans_annot)[grep("Patt",colnames(Pmeans_annot))]
  
  #Calculate correlation
  tmp <- cor(Pmeans_annot[,patterns], Pmeans_annot[,group_vector])
  
  #Melt correlation matrix   
  tmp_melt <- melt(tmp)
  colnames(tmp_melt) <- c("pattern", "group_vector","Correlation")
  tmp_melt$group_vector <- factor(tmp_melt$group_vector,levels=rev(group_vector[order(group_vector)]))
  
  #cluster patterns
  if(cluster_patterns== T) {
    if(is.null(cluster_patterns_vector)==F){
      tmp_melt$pattern = factor(tmp_melt$pattern,levels=unique(tmp_melt$pattern)[cluster_patterns_vector])
    } else {
      pattern.order<-order.dendrogram(as.dendrogram(pheatmap:::cluster_mat((tmp), distance=cluster_patterns_distance,method=cluster_patterns_method)))
      tmp_melt$pattern = factor(tmp_melt$pattern, levels=unique(tmp_melt$pattern)[pattern.order])
    }
  } else { tmp_melt$pattern = factor(tmp_melt$pattern, levels=rev(unique(tmp_melt$pattern)[order(tmp_melt$pattern,decreasing=T)])) }
  
  #cluster groups
  if(cluster_groups == T) {
    ref_table <- as.data.frame(group_vector)
    colnames(ref_table) <- ("group_vector")
    group.order<-order.dendrogram(as.dendrogram(pheatmap:::cluster_mat(t(tmp), distance=cluster_groups_distance,method=cluster_groups_method)))
    ref_table$group_vector <- factor(ref_table$group_vector, levels=ref_table$group_vector[group.order])
    #order_subgroups
    if(order_subgroups == T) {
    n_subgroups <- length(stringr::str_split(group_vector, pattern = " - ")[[1]])
    for(i in 1:n_subgroups) {
      ref_table[paste("subgroup",i,sep="")] <- as.factor(stringr::str_split_fixed(group_vector, pattern=" - ",
                                                                                  n=n_subgroups)[,i])
    }
    #order by subgroup1,2 ..
    for(i in n_subgroups:2) {
      ref_table <- ref_table[with(ref_table, order(ref_table[,i])),]
    }
    group_vector_ordered <- ref_table$group_vector
    tmp_melt$group_vector <-factor(tmp_melt$group_vector,levels=(group_vector_ordered)) 
    } else {
    tmp_melt$group_vector <-factor(tmp_melt$group_vector,levels=(ref_table$group_vector[group.order]))
    }
  } else { tmp_melt <- tmp_melt }
  
  #Check correlation range
  print(c(min(tmp_melt$`Correlation`),0, max(tmp_melt$`Correlation`)))
  
  #plot
  p<-ggplot(tmp_melt,aes(x=pattern,y=group_vector, fill=`Correlation`))
  p <- p + geom_tile(size=0.25) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
  p <- p + scale_fill_gradientn(colours =  c("dodgerblue","ghostwhite","brown2"),values=scales::rescale(c(min(tmp_melt$`Correlation`), 0,max(tmp_melt$`Correlation`))), limits=limits,  oob=scales::squish)
  if(square_tiles==T) {
    p <- p + coord_equal()
  } else { p <- p }
  return(p)
}