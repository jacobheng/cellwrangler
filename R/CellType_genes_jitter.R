#' Plot gene expression for all CellTypes as jitter plots
#'
#' @description CellType_genes_jitter is a wrapper function around monocle::plot_genes_jitter() that plots gene expression
#' as jitter plots facet_wrapped by CellType.
#' 
#' @param gene_to_plot official gene symbol of a gene to plot e.g. "Actb"
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @param cell_size the size (in points) of each cell used in the plot.
#' @param nrow the number of rows used when laying out the panels for each gene's expression.
#' @param ncol the number of columns used when laying out the panels for each gene's expression.
#' @param color the cell attribute (e.g. the column of pData(cds)) to be used to color each cell.
#' @param group the cell attribute (e.g. the column of pData(cds)) to group cells by on the horizontal axis.
#' @param plot_trend whether to plot a trendline tracking the average expression across the horizontal axis.
#' @param relative_expr Whether to transform expression into relative values.
#' @param color_scale a vector of colors to color barplots and error bars with.
#' @keywords gene_barplot()
#' @export
#' @return A ggplot2 object
#' @examples
#' gene_barplot("Actb", dat)
#' 
CellType_genes_jitter <- function(gene_to_plot ,cds, cell_size = 0.75, nrow=NULL, ncol=2, color="genotype", 
                                  group= "genotype", plot_trend=F, relative_expr = F, 
                                  color_scale = c("Black", "Red")){
  cds <- cds
  gene_id <- rownames(fData(cds)[fData(cds)$gene_short_name %in% c(gene_to_plot),])
  pData(cds)$CellType <- stringr::str_wrap(pData(cds)$CellType, width =5)
  plot_genes_jitter(cds[rownames(fData(cds)[gene_id,]),], cell_size = cell_size, nrow=nrow, ncol=ncol, 
                    color_by=color,grouping=group,relative_expr = relative_expr, plot_trend=plot_trend) + 
    facet_wrap("CellType") + ggtitle(gene_to_plot) + scale_color_manual(values = color_scale)
}