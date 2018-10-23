#'Normalize and log monocle cds
#'
#' log_normalize_cds() globally scales and log-transforms the expression matrix of a monocle CellDataSet object. 
#' The total count for each cell (represented by a column) is scaled to have the same total as all other cells 
#' (columns). A monocle cds with the log-normalized expression matrix is returned.
#'
#' @param cds a monocle CellDataSet object
#' @param filter_zeros logical - if genes with zero expression should be removed; defaults to TRUE
#' @param scale_to what value total read counts for each cell should be scaled to; by default, each cell's total
#' count is scaled to the median of all cells' totals; a numeric value can also be specified
#' @param logbase base for logarithm; defaults to 10; if NULL, no logarithm will be computed and the original
#' expression matrix will be returned
#' @keywords scale_log_cds
#' @export
#' @return a monocle CellDataSet object
#' @examples
#' dat_log <- log_normalize_cds(dat, scale_to = 10000, logbase=10)

log_normalize_cds <- function(cds, filter_zeros = TRUE, scale_to = "median", logbase=10) {
  if(filter_zeros == T) {
  nonzero_genes <- rownames(exprs(cds[rowSums(exprs(cds)) > 0,]))
  filtered_cds <- cds[nonzero_genes,]
  } else { filtered_cds <- cds }
  cell_totals <- Matrix::colSums(exprs(filtered_cds))
  if(scale_to == "median") {
  cell_median <- median(cell_totals)
  norm_exprs <- sweep(exprs(filtered_cds), 2, cell_median/cell_totals, "*")
  } else { if(is.numeric(scale_to) == T) {
    norm_exprs <- sweep(exprs(filtered_cds), 2, scale_to/cell_totals, "*")
  } else { message ("Need to specify scale factor!")}
  }
  if(logbase == NULL) {
    log_norm_exprs <- norm_exprs
  } else {
  log_norm_exprs <- log((norm_exprs+1), base = logbase)
  }
  new_cds <- newCellDataSet(as(as.matrix(log_norm_exprs), "sparseMatrix"),
                            phenoData=new("AnnotatedDataFrame",data=pData(filtered_cds)),
                            featureData=new("AnnotatedDataFrame",data=fData(filtered_cds)),
                            lowerDetectionLimit=1,
                            expressionFamily=negbinomial.size())
  return(new_cds)
}