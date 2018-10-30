#'Find expressed genes in each sample in CellDataSet
#'
#' find_sample_expressed_genes() splits an aggregated cds into separate samples based on the suffix at the 
#' end of each barcode. The function then proceeds to extract genes above the pre-specified UMI and cell 
#' expression threshold(s)
#' @param cds a monocle CellDataSet object
#' @param UMI_threshold UMI threshold for defining expressed genes. Genes with total UMIs equal to or above
#' this threshold will be selected.
#' @param cell_threshold cell threshold for defining expressed genes. Genes expressed in a total number of cells 
#' at or above this threshold will be selected.
#' @keywords find_sample_expressed_genes
#' @export
#' @return A list of vectors with each vector containing dispersed genes for each sample
#' @examples
#' find_sample_expressed_genes(dat,UMI_threshold = 1, cell_threshold = 10)


find_sample_expressed_genes <- function(cds, UMI_threshold = 1, cell_threshold = 1) {
  sample_numbers <- unique(str_split_fixed(colnames(exprs(cds)), 
                                           pattern = "-", n = 2)[, 2])
  expressed_genes_list <- lapply(sample_numbers, function(x, UMI_threshold, cell_thresholdl) {
    cds_subset <- cds[ , str_split_fixed(colnames(exprs(cds)), 
                                         pattern = "-", n = 2)[, 2] == x]
    expressed_genes_by_UMI <- Matrix::rowSums(exprs(cds_subset)) >= UMI_threshold
    cds_subset <- cds_subset[expressed_genes_by_UMI,]
    expressed_genes_by_cell <- Matrix::rowSums(exprs(cds_subset) > 0) >= cell_threshold
    cds_subset <- cds_subset[expressed_genes_by_cell,]
    expressed_genes <- rownames(cds_subset)
    return(expressed_genes)
  }, UMI_threshold = UMI_threshold, cell_threshold = cell_threshold)
  return(expressed_genes_list)
}