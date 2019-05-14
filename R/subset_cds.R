#' Subset cds by cell group
#'
#' @description subset_cds() subsets a CellDataSet object for a particular cell group and prepares the new CellDataSet
#' subset for the monocle differentialGeneTest() function. 
#' 
#' @param cell_group string specifying specific cell group e.g. "Rods"
#' @param cds a CellDataSet object used in the monocle package
#' @param group column name in pData(cds) specifying cell group that each cell belong to e.g. "CellType"
#' @param min_expr expression threshold to use for detectGenes() function
#' @keywords gene_barplot()
#' @export
#' @return CellDataSet object
#' @examples
#' subset_cds("Rods", dat, group = "CellType", min_expr = 1)


subset_cds <- function(cell_group, cds, group, min_expr=NULL) {
  cds_subset <- cds[, pData(cds)[,group] == cell_group]
  cds_subset <- estimateSizeFactors(cds_subset)
  cds_subset <- estimateDispersions(cds_subset)
  
  if(is.null(min_expr)) {
    min_expr <- tmp@lowerDetectionLimit
  } else { min_expr <- min_expr }
  cds_subset <- detectGenes(cds_subset, min_expr=min_expr)
  
  fData(cds_subset)$Prop_celltype_expressed <- Matrix::rowSums(exprs(tmp) > min_expr)/nrow(pData(tmp))
  fData(cds_subset)$Total_UMIs_per_gene_celltype <- Matrix::rowSums(exprs(cds_subset))
  fData(cds_subset)$Mean_celltype_UMIs <- Matrix::rowSums(exprs(cds_subset))/nrow(pData(cds_subset))
  fData(cds_subset)[,group] <- cell_group
  
  return(cds_subset)
}