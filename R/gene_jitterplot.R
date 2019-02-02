#' Plot expression of cells as jitter plot
#'
#' gene_jitterplot() plots the expression of a given gene as a jitter plot for a CellDataSet object.
#' 
#' @param genes a vector of gene name(s) to plot e.g. c("Actb", "Aldoa")
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @param group column name in pData(cds) to group cells by on the horizontal axis.
#' @param color column name in pData(cds) to color cells with.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param min_expr the minimum (untransformed) expression level to use in plotting the genes; expression 
#' values less than this threshold will be converted to 0
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param panel_order the order in which genes should be layed out (left-to-right, top-to-bottom)
#' @param plot_trend logical;whether to plot a trend line across groups in plot(s).
#' @param color_trend color to use for trend line.
#' @param label_by_short_name logical;label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param log_expr logical; whether to log-transform expression into log10(expression + 1)
#' @param relative_expr logical;Whether to transform expression into counts per 10,000
#' @keywords gene_jitterplot
#' @export
#' @return A ggplot2 object
#' @examples
#' gene_jitterplot(c("Actb", "Aldoa"), cds)

gene_jitterplot <- function (genes, cds, group= "genotype", color = NULL, cell_size = 0.75, min_expr = 0,
                             nrow = NULL, ncol = 1, panel_order = NULL, plot_trend = FALSE, 
                             color_trend = "orange", label_by_short_name = TRUE, 
                             log_expr= TRUE,relative_expr = FALSE) 

{
  cds_subset <- cds[cellwrangler::findGeneID(genes, cds), ]
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null( pData(cds_subset)$Total_UMI ) ) {
        stop("Error: to call this function with relative_expr=TRUE, you calculate Total_UMIs for each
             cell in pData")
      }
      cds_exprs <- (Matrix::t(Matrix::t(cds_exprs)/ pData(cds_subset)$Total_UMI )) * 10000
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr) == T) {
    min_expr <- 0
  } else { min_expr <- min_expr }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- 0
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                      levels = panel_order)
  }
  if(log_expr == T) {
    cds_exprs$log10_expression <- log10(cds_exprs$expression+1)
    p <- ggplot(aes_string(x = group, y = "log10_expression"), data = cds_exprs) + ylab( "Log10(Expression + 1) ")
  } else {
    p <- ggplot(aes_string(x = group, y = "expression"), data = cds_exprs) + ylab("Expression")
  }
 
  if (is.null(color) == FALSE) {
    p <- p + geom_jitter(aes_string(color = color), size = I(cell_size))
  }
  else {
    p <- p + geom_jitter(size = I(cell_size))
  }
  if (plot_trend == TRUE) {
    p <- p + stat_summary(color=color_trend, fun.data = "mean_cl_boot", 
                          size = 0.35)
    
    if(log_expr == T) { 
      p <- p + stat_summary(aes_string(x = group, y = "log10_expression", group = 1), color= color_trend, 
                            fun.data = "mean_cl_boot", size = 0.35, geom = "line") 
      } else {
      p <- p + stat_summary(aes_string(x = group, y = "expression", group = 1), color= color_trend, 
                          fun.data = "mean_cl_boot", size = 0.35, geom = "line") 
      }
  }
  p <- p + facet_wrap(~feature_label, nrow = nrow, ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    p <- p + expand_limits(y = c(min_expr, 1))
  }
  p <- p + xlab(group)
  p <- p + monocle:::monocle_theme_opts()
  p
}