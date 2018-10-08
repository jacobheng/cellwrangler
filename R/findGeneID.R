#' Find Ensembl gene IDs given gene names
#'
#' Using gene names as input, findGeneID() returns the corresponding Ensembl gene IDs in the fData
#' of the specified CellDataSet object.
#' @param gene_names a vector of official gene names e.g. "Acta2"
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @keywords findGeneID
#' @export
#' @examples
#' findGeneID(c("Acta2","Cdh5","Lum"), myCDS)


findGeneID<-function(gene_names, cds){
  GeneIDs <- rownames(fData(cds))[fData(cds)$gene_short_name %in% gene_names]
  GeneIDs <- unique(GeneIDs)
  return(GeneIDs)
}