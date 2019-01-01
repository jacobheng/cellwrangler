#'Determine if cells are above an expression threshold of a gene
#'
#' dim_coord_barcode() returns the rownames of cells based on the coordiates of a dimensionality reduction
#' projection (e.g. t-SNE, UMAP).
#' 
#' @param cds a monocle CellDataSet object.
#' @param gene_name official gene name e.g. "Actb".
#' @param threshold_value numeric; threshold of expression of gene.
#' @keywords exprs_threshold
#' @export
#' @return A boolean vector
#' @examples
#' exprs_threshold(dat, "Actb", 1)

#exprs_threshold function
exprs_threshold <- function(cds,gene_name,threshold_value) {
  tmp <- exprs(cds)[findGeneID(gene_name, cds),] >= threshold_value
  return(tmp)
}
