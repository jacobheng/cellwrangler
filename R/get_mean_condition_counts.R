#' Get mean read counts for cells belonging to a condition
#'
#' @description get_mean_condition_counts finds the mean counts for cells belonging to a condition e.g. genotype 
#' in a CellDataSet object and returns the values in the pData of CellDataSet object
#' 
#' @param cds a CellDataSet object used in the monocle package
#' @param condition column name in pData(cds) that specifies the condition for classifying cells e.g. genotype
#' @param read_type string specifying type of reads in expression matrix; defaults to "UMIs"
#' @keywords get_mean_condition_counts
#' @export
#' @return CellDataSet object
#' @examples
#' dat <- get_mean_condition_counts(dat, condition="genotype")

get_mean_condition_counts <- function(cds, condition, read_type = "UMIs") {
  tmp <- cds
  conditions <- unique(pData(cds)[,condition])
  
  for(i in 1:length(conditions)) {
    fData(tmp)[paste0("Mean_", conditions[i], "_", read_type)] <- Matrix::rowSums(exprs(tmp[,pData(tmp)[,condition]==conditions[i]])) / sum(pData(tmp)[,condition] == conditions[i])
  }
  return(tmp)
 
}