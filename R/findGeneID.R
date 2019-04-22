#' Find Ensembl gene IDs given gene names
#'
#' @description Using gene names as input, findGeneID() returns the corresponding Ensembl gene IDs in the fData
#' of the specified CellDataSet object.
#' @param gene_names a vector of official gene names e.g. "Acta2"
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @keywords findGeneID
#' @export
#' @examples
#' findGeneID(c("Acta2","Cdh5","Lum"), myCDS)


findGeneID<-function(gene_names, cds){
  
  gene_ref <- as.data.frame(as.character(gene_names))
  colnames(gene_ref) <- c("gene_short_name")
  gene_ref$gene_id <- rownames(fData(cds))[match(gene_ref$gene_short_name, fData(cds)$gene_short_name)]
  
  GeneIDs <- as.character(gene_ref$gene_id)
  
  return(GeneIDs)
}
