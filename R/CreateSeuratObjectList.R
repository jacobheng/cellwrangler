#'Create a list of Seurat Objects
#'
#' CreateSeuratObjectList() creates a list of Seurat objects for use in the RunMultiCCA function in the R package
#' Seurat. The function splits a single aggregated barcoded expression matrix into separate matries based on the 
#' suffix at the end of each cellular barcode e.g. ("ACTAGGAGACAGGTGC-1").
#'  
#' @param exprs_matrix an expression matrix; row names must be gene names and column names must be cell barcodes
#' @param phenoData dataframe where rows are cells and rownames are cell barcodes, columns are fields defining 
#' cellular characteristics. Row names of phenoData must match column names of exprs_matrix. 
#' @keywords CreateSeuratObjectList
#' @export
#' @return A list of Seurat Objects
#' @examples
#' list <- CreateSeuratObjectList(exprs_matrix = exprs(dat), phenoData = pData(dat))


CreateSeuratObjectList <- function(exprs_matrix, phenoData) {
  sample_numbers <- unique(str_split_fixed(colnames(exprs_matrix), 
                                           pattern = "-", n = 2)[, 2])
  if(colnames(exprs_matrix) == rownames(phenoData)) {
    message("colnames of expression matrix and rownames phenoData match")
  } else {message("Error: colnames of expression matrix and rownames phenoData do not match!")}
  
  SeuratObjectList <- lapply(sample_numbers, function(x) {
    exprs_subset <- exprs_matrix[ , str_split_fixed(colnames(exprs_matrix), 
                                         pattern = "-", n = 2)[, 2] == x]
    phenoData_subset <- phenoData[str_split_fixed(rownames(phenoData), 
                                                  pattern = "-", n = 2)[, 2] == x,]
    SeuratObj <- CreateSeuratObject(raw.data = exprs_subset, meta.data=phenoData_subset)
    SeuratObj <- NormalizeData(SeuratObj)
    SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
    SeuratObj <- ScaleData(SeuratObj)
    
    return(SeuratObj)
  })
  return(SeuratObjectList)
}