#'Split a cellranger-aggregated matrix into multiple matrices based on barcodes
#'
#' @description split_aggr_matrix() splits an expression matrix derived from a directory of 10X Genomics scRNAseq 
#' data aggregated from multiple samples created using the cellranger pipeline. Based on the unique numeric 
#' value assigned to barcodes of each sample by the cellranger aggr function.
#'
#' @param aggr_mtx an expression matrix derived from the cellranger aggr function
#' @keywords split_aggr_matrix
#' @export
#' @return a dgTMatrix 
#' @examples
#' sample_raw_mtx <- split_aggr_mtx(raw_mtx)

split_aggr_mtx <- function(aggr_mtx) {
  sample_numbers <- unique(str_split_fixed(colnames(aggr_mtx),pattern="-",n=2)[,2])
  aggr_mtx_split <- lapply(sample_numbers, function(x){
    tmp <- aggr_mtx[ ,str_split_fixed(colnames(aggr_mtx),pattern="-",n=2)[,2]==x ]
    return(tmp)
  } )
  return(aggr_mtx_split)
}