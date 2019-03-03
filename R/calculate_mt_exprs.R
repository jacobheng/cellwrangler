#'Calculate mitochondrial gene expression
#'
#' calculate_mt_RNAs() calculates the expression of mitochondrial genes for a monocle CellDataSet
#' object.
#' @param cds monocle CellDataSet object
#' @param genome genome to which transcripts are aligned to. Accepts either "mouse" or "human".
#' @param return_prop logical; whether to return the value of mitochondrial transcripts as 
#' a proportion of total transcripts for each cell; resulting output will be a dataframe. 
#' If FALSE, will return only absolute number of transcripts per cell.
#' @keywords mt_genes, mt_RNAs
#' @export
#' @return A vector or a dataframe
#' @examples
#' calculate_mt_exprs




calculate_mt_exprs <- function(cds, genome = "mouse", return_prop = F) {

  if(genome == "mouse") {
  mt_genes<- base::grepl("^mt-",fData(cds)$gene_short_name)
  } else {
    if(genome == "human") {
  mt_genes<- base::grepl("^MT-",fData(cds)$gene_short_name)
    } else { warning("Need to specify genome as human or mouse!") }
  }
  if(return == F) {
    mt_RNAs <- Matrix::colSums(exprs(cds[mt_genes,])) 
    tmp <- mt_RNAs
    names(tmp) <- rownames(pData(cds))
  } else {
    mt_RNAs <- Matrix::colSums(exprs(cds[mt_genes,])) 
    Total_RNAs <- Matrix::colSums(exprs(cds))
    mt_prop <- mt_RNAs/Total_RNAs
    tmp <- as.data.frame(cbind(mt_RNAs, mt_prop))
    rownames(tmp) <- rownames(pData(cds))
  }
  
  return(tmp) 
}



