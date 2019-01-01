#'Plot a set of patterns on a dimensionality reduction plot
#'
#' plot_patternSet_dim() plots Pmeans of a CoGAPS pattern set on a dimensonality reduction plots.
#'
#' @param CoGAPS_res_set a CoGAPS result set. If set to NULL, a Pmeans matrix may be supplied directly as input
#' for the Pattern_set parameter.
#' @param Pattern_set a string specifying the name of a pattern set (e.g. "nP20") in a CoGAPS_res_set. 
#' Alternatively, if CoGAPS_res_set is set to NULL, a Pmeans matrix may be supplied directly.
#' @param annotation an annotation dataframe. Must contain columns with column names corresponding to values
#' in the group_vector. 
#' @param dim_reduction a dataframe specifying the coordinates of a dimensionality reduction method (e.g. 
#' t-SNE, UMAP etc.). Each column should specify one dimension. Rownames should correspond to the rownames
#' of the annotation object. Colnames will be used as axis labels.
#' @param facet_wrap_by a string specifying variable for facet_wrap; defaults to NULL.
#' 
#' @keywords plot_pattern_dim
#' @export
#' @return a list of ggplot objects
#' @examples
#' plot_pattern_dim(CoGAPS_res_set= my_CoGAPS_res, Pattern_set="nP30",  annotation=pData(cds), 
#' dim_reduction = UMAP_proj)





#plot_patterns_dim function
plot_patternSet_dim <- function(Pattern_set, CoGAPS_res, annotation, dim_reduction, facet_wrap_by=NULL) {
  if(is.null(CoGAPS_res_set)==F) {
    Pmeans <- (t(CoGAPS_res_set[[Pattern_set]]$Pmean))} else { Pmeans <- Pattern_set }
  
  Pmeans <- Pmeans[unique(rownames(Pmeans)),]
  Pmeans_annot <- merge_by_rownames(annotation, Pmeans, all.x = F, all.y=T)
  #Get dim_reduction names
  dim_reduction_names <- colnames(dim_reduction)
  colnames(dim_reduction) <- c("Dimension1","Dimension2")
  #Merge dim reduction with Pmeans_annot
  Pmeans_annot <- merge_by_rownames(Pmeans_annot , dim_reduction, all.x = T, all.y=F)
  for(i in c( (ncol(annotation)+1):(ncol(Pmeans)+ncol(annotation)) )){
    Pmeans_annot[,i]<-as.numeric(Pmeans_annot[,i])
  }
  

  #Plot by dim_reduction
  Pmeans_plots <- lapply(c(1:ncol(Pmeans)), function(x){
    pattern_to_plot <- ncol(annotation) + x
    p <- ggplot(Pmeans_annot, aes(x=Dimension1,y=Dimension2)) + 
      geom_point(color=Pmeans_annot[,pattern_to_plot]) +
      scale_colour_gradientn(colours =c("grey90","darkred"),
                             values=scales::rescale(c(min(Pmeans_annot[,pattern_to_plot]),
                                                      max(Pmeans_annot[,pattern_to_plot])))) +
      xlab(dim_reduction_names[1]) + ylab(dim_reduction_names[2]) + 
      guides(fill=guide_legend(title=paste("Pattern",pattern_num,sep=""))) + theme_classic
    
    if(is_null(facet_wrap_by) == F) {
      p <- p + facet_wrap(facet_wrap_by)
    } else { p <- p }
    return(p)
    })
  
  #Get all plots
  return(Pmeans_plots)
}