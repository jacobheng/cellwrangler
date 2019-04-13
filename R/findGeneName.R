#' Find gene names given Ensembl gene IDs
#'
#' @description Using Ensembl gene IDs as input, findGeneName() returns the corresponding gene names in the fData
#' of the specified CellDataSet object.
#' @param gene_IDs a vector of Ensembl gene IDs e.g ENSMUSG00000035783
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @keywords findGeneName
#' @export
#' @return string(s)
#' @examples
#' findGeneName(c("ENSMUSG00000035783","ENSMUSG00000031871","ENSMUSG00000036446"), myCDS)


findGeneName<-function(gene_IDs, cds){
  
  gene_ref <- as.data.frame(as.character(gene_IDs))
  colnames(gene_ref) <- c("gene_id")
  gene_ref$gene_short_name <- fData(cds)$gene_short_name[match(gene_ref$gene_id, rownames(fData(cds)))]
  
  GeneNames <- as.character(gene_ref$gene_short_name)
  
  return(GeneNames)
}