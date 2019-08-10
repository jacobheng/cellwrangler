#' Plot mean gene expression as bar plot
#'
#' @description gene_barplot() plots the mean expression of a given gene as a bar plot with standard error bars.
#' 
#' @param genes a vector of gene name(s) to plot e.g. c("Actb", "Aldoa")
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @param group column name in pData(cds) to group the bar plots by on the horizontal axis.
#' @param color column name in pData(cds) to color the bar plots with.
#' @param facet_wrap a column in pData(cds) to facet_wrap the plots with.
#' @param facet_genes only used if facet_wrap is NULL and there are multiple genes being plotted; if set to "cols", will
#' facet genes into columns; otherwise will facet genes into rows
#' @param plot_trend logical;whether to plot a trend line across groups in plot(s).
#' @param color_trend color to use for trend line.
#' @keywords gene_barplot()
#' @export
#' @return A ggplot2 object
#' @examples
#' gene_barplot(c("Actb", "Aldoa"), dat)

gene_barplot <- function (genes, cds, group = "genotype", color = "genotype", facet_wrap = NULL, facet_genes="cols", 
                          plot_trend = F, color_trend = "orange") 
{
  cds_subset <- cds[cellwrangler::findGeneID(genes, cds), ]
  exprs_values <- Matrix::t(as.matrix(exprs(cds_subset)))
  colnames(exprs_values) <- cellwrangler::findGeneName(colnames(exprs_values), cds)  
  tmp <- merge_by_rownames(pData(cds), exprs_values)
  
  if(length(genes) == 1) {
    
    colnames(tmp) <- c(colnames(pData(cds)), "exprs")
    p <- ggplot(tmp, aes_string(x = group, y = "exprs")) + ggtitle(genes)
    if(is.null(facet_wrap)==FALSE) {
      p <- p + facet_wrap(facet_wrap)
    } else { p <- p }
    
  } else {
    
    tmp_melt <- melt(tmp, id.vars= colnames(pData(cds)))
    colnames(tmp_melt) <- c(colnames(pData(cds)), "gene","exprs")
    print(head(tmp_melt))
    p <- ggplot(tmp_melt, aes_string(x = group, y = "exprs")) 
    if(is.null(facet_wrap)==FALSE) {
      p <- p + facet_grid(paste("gene~",facet_wrap,sep=""))
    } else { 
      if(facet_genes == "cols") {
       p <- p + facet_grid(cols = vars(gene))  
      } else { p <- p + facet_grid(rows = vars(gene)) }
      }
    
  }
  p <- p  + stat_summary(fun.y = mean, geom = "bar", aes_string(fill = color)) + 
    stat_summary(fun.data = mean_se, geom = "errorbar", aes_string(color = color), width = 0.25) + 
    ylab("Mean UMIs per cell")
  
  if (plot_trend == T) {
    p <- p + stat_summary(color=color_trend, fun.data = "mean_cl_boot", 
                          size = 0.35)
    p <- p + stat_summary(aes_string(x = group, y = "exprs", group = 1), color=color_trend, fun.data = "mean_cl_boot", 
                          size = 0.35, geom = "line")
  }
  else {
    p <- p
  }
  
  return(p)
}