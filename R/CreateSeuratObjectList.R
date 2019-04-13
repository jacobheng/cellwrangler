#'Create a list of Seurat Objects
#'
#' @description CreateSeuratObjectList() creates a list of Seurat objects for use in the RunMultiCCA function in the R package
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


CreateSeuratObjectList <- function(exprs_matrix, phenoData, scaleData = FALSE) {
  sample_numbers <- unique(str_split_fixed(colnames(exprs_matrix), 
                                           pattern = "-", n = 2)[, 2])
  
  SeuratObjectList <- lapply(sample_numbers, function(x, scale) {
    exprs_subset <- exprs_matrix[ , str_split_fixed(colnames(exprs_matrix), 
                                         pattern = "-", n = 2)[, 2] == x]
    exprs_subset <- as(as.matrix(exprs_subset), "sparseMatrix")
    phenoData_subset <- phenoData[str_split_fixed(rownames(phenoData), 
                                                  pattern = "-", n = 2)[, 2] == x,]
    SeuratObj <- CreateSeuratObject(raw.data = exprs_subset, meta.data=phenoData_subset)
    SeuratObj <- NormalizeData(SeuratObj)
    SeuratObj <- FindVariableGenes(SeuratObj, do.plot = F, display.progress = F)
    if(scale == T) {
    SeuratObj <- ScaleData(SeuratObj)
    } else { SeuratObj <- SeuratObj }
    return(SeuratObj)
  }, scale = scaleData)
  return(SeuratObjectList)
}