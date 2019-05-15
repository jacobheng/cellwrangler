#' Find genes expressed in a CellDataSet object
#'
#' @description get_expressed_genes() finds genes expressed in a CellDataSet object according to defined
#' parameters.
#' 
#' @param cds a CellDataSet object used in the monocle package
#' @param cell_threshold numerical threshold for number of cells that must express a particular gene e.g. 1; 
#' if min_prop is specified, this number will be overridden. 
#' @param min_prop minimum proportion of cells that must express a particular gene e.g. 0.05 for the gene
#' to be considered expressed
#' @param min_condition a column in pData(cds) that specifies a condition e.g. genotype to which the min_prop 
#' argument will be applied; must be specified together with the min_prop argument. Takes the minimum proportion
#' of the condition with fewer cells as the threshold for finding expressed genes.
#' @keywords get_expressed_genes
#' @export
#' @return CellDataSet object
#' @examples
#' expressed_genes <- get_expressed_genes(dat, min_prop=0.05, min_condition="genotype")


get_expressed_genes <- function(cds, cell_threshold=1, min_prop=NULL, min_condition=NULL) {
  cds <- cds
  if(is.null(min_prop) == F) {
    if(is.null(min_condition) == F) {
      numCellThreshold<- min_prop*min(table(pData(cds)[,condition]))
    } else {
      numCellThreshold <- min_prop*nrow(pData(cds))
    }
  } else { numCellThreshold <- cell_threshold }
  
  tmp <- row.names(subset(fData(cds), num_cells_expressed >= numCellThreshold))
  return(tmp)
}
