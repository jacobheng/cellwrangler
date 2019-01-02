#' Plot biological coefficent of variation and annotate genes
#'
#' plot_bcv() plots the mean - variance (bcv) relationship, as well as colors and annotates genes of interest.
#' 
#' @param annotated_fData annotated fData object produced by the cellwrangler annotate_bcv() function.
#' @param genes_to_color vector of Ensembl gene ids to color
#' @param color_scale color scale used for coloring genes specified in genes_to_color.
#' @param genes_to_annotate vector of Ensembl gene ids to annotate
#' @keywords plot_bcv
#' @export
#' @examples
#' plot_bcv(annotated_fData)



plot_bcv <- function(annotated_fData, genes_to_color, color_scale=c('grey50','red'), genes_to_annotate, 
                     title = NULL) {
  
  annotated_fData$to_color <-FALSE
  annotated_fData$to_color[annotated_fData$id %in% genes_to_color]<-TRUE
  
  p <- ggplot(data=annotated_fData)
  p <- p + geom_point(aes(x=log(mean_exprs),y=as.vector(log(bcv)),color=to_color),alpha=0.5,size=0.6) + 
    scale_color_manual(values= color_scale ) +
    geom_smooth(aes(x=log(mean_exprs),y=as.vector(log(bcv))),color="black",method="auto",se=TRUE) + 
    geom_text_repel(aes(x=log(mean_exprs),y=as.vector(log(bcv)),label=gene_short_name),color="black",size=3,
                    data=subset(annotated_fData,id %in% genes_to_annotate)) +
    ggtitle(title) + guides(color=FALSE) + ylab("log(bcv)") + theme_bw()
  
  return(p)
  
}