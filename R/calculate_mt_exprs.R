#'Calculate mitochondrial gene expression
#'
#' @description calculate_mt_exprs() calculates the expression of mitochondrial genes for a monocle CellDataSet
#' object.
#' @param cds monocle CellDataSet object
#' @param genome genome to which transcripts are aligned to. Accepts either "mouse" or "human".
#' @param return_prop logical; whether to return the value of mitochondrial transcripts as 
#' a proportion of total transcripts for each cell; resulting output will be the input CDS with
#' the values of mt_RNAs and mt_prop attached to the corresponding pData. If FALSE, will only
#' return the absolute number of transcripts per cell as a vector. Defaults to TRUE.
#' @keywords mt_genes, mt_RNAs
#' @export
#' @return A vector or a cds (see above)
#' @examples
#' calculate_mt_exprs




calculate_mt_exprs <- function(cds, genome = "mouse", return_prop = T) {

  if(genome == "mouse") {
  mt_genes<- base::grepl("^mt-",fData(cds)$gene_short_name)
  } else {
    if(genome == "human") {
  mt_genes<- base::grepl("^MT-",fData(cds)$gene_short_name)
    } else { warning("Need to specify genome as human or mouse!") }
  }
  if(return_prop == F) {
    mt_RNAs <- Matrix::colSums(exprs(cds[mt_genes,])) 
    tmp <- mt_RNAs
    names(tmp) <- rownames(pData(cds))
  } else {
    mt_RNAs <- Matrix::colSums(exprs(cds[mt_genes,])) 
    Total_RNAs <- Matrix::colSums(exprs(cds))
    mt_prop <- mt_RNAs/Total_RNAs
    tmp <- cds
    pData(tmp)$mt_RNAs <- mt_RNAs
    pData(tmp)$mt_prop <- mt_prop
  }
  
  return(tmp) 
}



