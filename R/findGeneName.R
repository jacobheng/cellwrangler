#' Find gene names given Ensembl gene IDs
#'
#' Using Ensembl gene IDs as input, findGeneName() returns the corresponding gene names in the fData
#' of the specified CellDataSet object.
#' @param gene_IDs a vector of Ensembl gene IDs e.g ENSMUSG00000035783
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @param unique logical; whether to return only unique gene names
#' @keywords findGeneName
#' @export
#' @return string(s)
#' @examples
#' findGeneName(c("ENSMUSG00000035783","ENSMUSG00000031871","ENSMUSG00000036446"), myCDS)


findGeneName<-function(gene_IDs, cds, unique = T){
  GeneNames <- fData(cds[rownames(cds) %in% gene_IDs,])$gene_short_name
  if(unique == T) {
  GeneNames <- unique(GeneNames) } else {
  GeneNames <- GeneNames
  }
  return(GeneNames)
}