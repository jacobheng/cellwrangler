#'Find dispersed genes in each sample in CellDataSet
#'
#' @description find_sample_dispersed_genes() splits an aggregated cds into separate samples based on the suffix 
#' at the end of each barcode and estimates SizeFactors and Dispersions for each sample using the package monocle. 
#' The function then proceeds to extract genes above pre-specified mean expression and dispersion empirical
#' thresholds
#' @param cds a monocle CellDataSet object
#' @param mean_exprs mean expression threshold for defining dispersed genes
#' @param disp_empirical dispersion empirical threshold for defining dispersed genes
#' @keywords find_sample_dispersed_genes
#' @export
#' @return A list of vectors with each vector containing dispersed genes for each sample
#' @examples
#' find_sample_dispersed_genes(dat, mean_exprs=0, disp_empirical=1)


find_sample_dispersed_genes <- function(cds, mean_exprs = 0, disp_empirical=1) {
  sample_numbers <- unique(str_split_fixed(colnames(exprs(cds)), 
                                           pattern = "-", n = 2)[, 2])
  cds_split <- lapply(sample_numbers, function(x, mean_exprs, disp_empirical) {
    cds_subset <- cds[ , str_split_fixed(colnames(exprs(cds)), 
                                         pattern = "-", n = 2)[, 2] == x]
    cds_subset <- estimateSizeFactors(cds_subset)
    cds_subset <- estimateDispersions(cds_subset)
    disp_table <- dispersionTable(cds_subset)
    dispersed_genes <- subset(disp_table, mean_expression >= mean_exprs & dispersion_empirical >= disp_empirical*dispersion_fit)$gene_id
    return(dispersed_genes)
  }, mean_exprs = mean_exprs, disp_empirical=disp_empirical)
  return(cds_split)
}