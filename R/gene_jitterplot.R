#' Plot expression of cells as jitter plot
#'
#' gene_jitterplot() plots the expression of a given gene as a jitter plot for a CellDataSet object.
#' 
#' @param genes a vector of gene name(s) to plot e.g. c("Actb", "Aldoa")
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @param group column name in pData(cds) to group cells by on the horizontal axis.
#' @param color column name in pData(cds) to color cells with.
#' @param cell_size the size (in points) of each cell used in the plot
#' @param min_expr the minimum (untransformed) expression level to use in plotted the genes.
#' @param nrow the number of rows used when laying out the panels for each gene's expression
#' @param ncol the number of columns used when laying out the panels for each gene's expression
#' @param plot_trend logical;whether to plot a trend line across groups in plot(s).
#' @param color_trend color to use for trend line.
#' @param label_by_short_name logical;label figure panels by gene_short_name (TRUE) or feature id (FALSE)
#' @param relative_expr logical;Whether to transform expression into relative values
#' @keywords gene_jitterplot
#' @export
#' @return A ggplot2 object
#' @examples
#' gene_jitterplot(c("Actb", "Aldoa"), cds)



gene_jitterplot <- function (genes, cds, group= "genotype", color = NULL, cell_size = 0.75, min_expr = NULL,
                             nrow = NULL, ncol = 1, panel_order = NULL, plot_trend = FALSE, 
                             color_trend = "orange", label_by_short_name = TRUE, relative_expr = TRUE) 

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
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
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
  p <- ggplot(aes_string(x = group, y = "expression"), data = cds_exprs)
  if (is.null(color) == FALSE) {
    p <- p + geom_jitter(aes_string(color = color), size = I(cell_size))
  }
  else {
    p <- p + geom_jitter(size = I(cell_size))
  }
  if (plot_trend == TRUE) {
    p <- p + stat_summary(color=color_trend, fun.data = "mean_cl_boot", 
                          size = 0.35)
    p <- p + stat_summary(aes_string(x = group, y = "expression", group = 1), color= color_trend, 
                          fun.data = "mean_cl_boot", size = 0.35, geom = "line")
  }
  p <- p + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    p <- p + expand_limits(y = c(min_expr, 1))
  }
  p <- p + ylab("Expression") + xlab(group)
  p <- p + monocle:::monocle_theme_opts()
  p
}