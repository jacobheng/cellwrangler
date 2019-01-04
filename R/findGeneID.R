#' Find Ensembl gene IDs given gene names
#'
#' Using gene names as input, findGeneID() returns the corresponding Ensembl gene IDs in the fData
#' of the specified CellDataSet object.
#' @param gene_names a vector of official gene names e.g. "Acta2"
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @param unique logical; whether to return only unique gene ids
#' @keywords findGeneID
#' @export
#' @examples
#' findGeneID(c("Acta2","Cdh5","Lum"), myCDS)


findGeneID<-function(gene_names, cds, unique = T){
  GeneIDs <- rownames(fData(cds)[fData(cds)$gene_short_name %in% gene_names,])
  if(unique == T) {
  GeneIDs <- unique(GeneIDs) } else {
  GeneIDs <- GeneIDs
  }
  return(GeneIDs)
}
