#' Plot mean gene expression as bar plot
#'
#' gene_barplot() plots the mean expression of a given gene as a bar plot with standard error bars.
#' 
#' @param gene_to_plot official gene symbol of a gene to plot e.g. "Actb"
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @param group column name in pData(cds) to group the bar plots by on the horizontal axis.
#' @param color column name in pData(cds) to color the bar plots with.
#' @param color_scale a vector of colors to color barplots and error bars with.
#' @keywords gene_barplot()
#' @export
#' @return A ggplot2 object
#' @examples
#' gene_barplot("Actb", dat)

gene_barplot <- function(gene_to_plot, cds, group = "genotype", color = "genotype", 
                         color_scale = c("Red", "Black")){
  tmp <- merge_by_rownames((exprs(cds)[rownames(fData(cds)[fData(cds)$gene_short_name %in% gene_to_plot,]),]), 
               pData(cds))
  colnames(tmp) <- c("exprs", colnames(pData(cds)))
  #plot
  p <- ggplot(tmp, aes_string(x = group, y = "exprs")) +
    stat_summary(fun.y = mean, geom = "bar", aes_string(fill=color))  +
    stat_summary(fun.data = mean_se, geom = "errorbar", aes_string(color=color), width = 0.25) + 
    ylab("Mean UMIs per cell") + scale_fill_manual(values = color_scale) + 
    scale_color_manual(values = color_scale) +
    ggtitle(gene_to_plot)
  
  return(p)
}