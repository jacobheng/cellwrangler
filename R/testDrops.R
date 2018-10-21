#'Determine whether droplets are cells or blank drops in an expression matrix
#'
#' testDrops() is a wrapper function around the emptyDrops() function in the excellent SoupX package. Returns
#' a list object consisting of a results dataframe from the emptyDrops() function, as well as barcodes 
#' of cells and empty droplets.
#'
#' @param raw_exprs_mtx a raw expression matrix created with the cellranger pipeline; rows are genes and 
#' columns are barcodes
#' @param lower_UMI_threshold numeric; parameter for the testDrops function. All barcodes with total UMI counts
#' below this value are considered empty droplets (i.e. not cells)
#' @param upper_UMI_threshold numeric; parameter for the testDrops function. All barcodes with total UMI counts
#' above this value are considered cells
#' @param test.ambient A logical scalar indicating whether results should be returned for barcodes with totals 
#' less than or equal to lower_UMI_threshold.
#' @param FDR numeric; parameter for the testDrops function. FDR value to use in the testDrops function when
#' statistically determining cells vs non-cells
#' @param niters An integer scalar specifying the number of iterations to use for the Monte Carlo p-value 
#' calculations.
#' @keywords testDrops
#' @export
#' @return a testDrops object with a results dataframe, cell_barcodes, and blank_barcodes
#' @examples
#' testDrops_res <- testDrops(raw_exprs_mtx, lower_UMI_threshold=100, upper_UMI_threshold=500, 
#' test.ambient=T, FDR=0.01, niters=10000)

testDrops <- function(raw_exprs_mtx, lower_UMI_threshold=100, upper_UMI_threshold=500, 
                      test.ambient=T, FDR=0.01, niters=10000) {
  set.seed(lower_UMI_threshold)
  tmp <- DropletUtils::emptyDrops(raw_exprs_mtx, lower=lower_UMI_threshold, retain=upper_UMI_threshold,
                                  test.ambient=test.ambient, niters=niters)
  tmp$is_cell <- tmp$FDR <= FDR
  tmp <- as.data.frame(tmp)
  cell_barcodes <- rownames(tmp[tmp$is_cell == TRUE & tmp$Total > lower_UMI_threshold,])
  blank_barcodes <- rownames(tmp[tmp$is_cell == FALSE & tmp$Total < upper_UMI_threshold,])
  res <- list(data=tmp, cell_barcodes=cell_barcodes, blank_barcodes=blank_barcodes)
  return(res)
}