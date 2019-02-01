#' Plot gene expression in cells
#'
#' plot_gene_exprs() plots levels of expression of a specified gene or genes as a color gradient. A 
#' two-dimensional dimensionality reduction object with coordinates for each cell (e.g. by t-SNE or UMAP dimensionality 
#' reduction) is required. 
#' @param cds a CellDataSet object or equivalent
#' @param genes a vector of gene name(s) to plot
#' @param dim_reduction a dataframe specifying the coordinates of a dimensionality reduction output (e.g. 
#' t-SNE, UMAP etc.). Each column should specify one dimension.
#' @param color_scale a vector of two colors for color gradient
#' @param limits a vector containing two number specifying the lower and upper limits of the color scale 
#' respectively.
#' @param rescale logical; whether to rescale color scale as a continuous gradient from the minimum to 
#' the maximum value of expression for specified gene; defaults to FALSE
#' @param cell_size, the size of the point representing each cell
#' @param group_vector a vector classifying cells e.g. genotype or condition
#' @param group_facet logical; whether to facet expression plots by cell groups specified in the group
#' parameter. Defaults to FALSE. If TRUE, a group vector must be supplied
#' @param title title of plot
#' @keywords plot_gene_exprs
#' @export
#' @return A ggplot2 object
#' @examples
#' plot_gene_exprs(myCDS, genes=c("Acta2","Myl9"), dim_reduction=tsne_proj[,c("TSNE.1", "TSNE.2")],
#' color_scale=c("blue","red"), limits = c(0,5), rescale= TRUE, cell_size =0.1, group_vector=pData(dat)$genotype,
#' group_facet=TRUE, title=NULL )

plot_gene_exprs <- function (cds, genes, dim_reduction, color_scale=c("slategray1","red"), limits = c(0, 10), 
                             rescale=F, cell_size = 0.1,  group_vector= NULL, group_facet=NULL, 
                                  title = NULL) 
{
  cds_subset <- cds[cellwrangler::findGeneID(genes,cds),]
  exprs_values <- t(as.matrix(exprs(cds_subset)))
  colnames(exprs_values) <- cellwrangler::findGeneName(colnames(exprs_values),cds)
  dim_reduction_names <- colnames(dim_reduction)
  colnames(dim_reduction) <- c("Component.1", "Component.2")
  if(is.null(group_vector) == F) {
    group_vector <- as.data.frame(group_vector)
    colnames(group_vector) <- c("group")
    proj_exprs <- data.frame(cbind(dim_reduction, group_vector, exprs_values))
    proj_exprs_melt <- melt(proj_exprs, id.vars = c("Component.1", "Component.2","group"))
    print(head(proj_exprs_melt))
  } else {
    proj_exprs <- data.frame(cbind(dim_reduction, exprs_values))
    proj_exprs_melt <- melt(proj_exprs, id.vars = c("Component.1", "Component.2"))
    print(head(proj_exprs_melt))
  }
  if(!is.null(group_facet)){
    p <- ggplot(proj_exprs_melt, aes(Component.1, Component.2)) + 
      geom_point(aes(colour = value), size = cell_size) + 
      facet_grid(group_facet)
  } else {
    p <- ggplot(proj_exprs_melt, aes(Component.1, Component.2)) + 
      geom_point(aes(colour = value), size = cell_size) + 
      facet_wrap(~variable) }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if(rescale==T) {
    p <- p + scale_color_gradientn(name= "Expression" ,colours =  color_scale, 
                                   values=scales::rescale(c(min(proj_exprs_melt$value),
                                                            max(proj_exprs_melt$value))))
  } else {
    p <- p + scale_color_gradient(name= "Expression",low=color_scale[1],high=color_scale[2],
                                  limits=limits,oob=scales::squish)
  }
  p <- p + theme_bw() + xlab(dim_reduction_names[1]) + ylab(dim_reduction_names[2]) + theme(plot.title = element_text(hjust = 0.5),
                                                                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                      strip.background=element_blank())
  return(p)
}