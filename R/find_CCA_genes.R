#'Find high variance genes for CCA
#'
#' @description find_CCA_genes() defines genes to use for CCA by finding high-variance genes in a list of Seurat objects
#' that are highly variable in a minimum number of samples/objects.
#'  
#' @param SeuratObjectList a list of Seurat objects
#' @param num_hvg_genes number of high-variance genes per sample to examine starting from the most variable
#' @param num_samples_variable minimum numher of samples the high-variance genes must be highly variable in
#' @keywords CreateSeuratObjectList
#' @export
#' @return A list of Seurat Objects
#' @examples
#' genes.use <- find_CCA_genes(SeuratObList,  num_hvg_genes=1000, num_samples_variable=2)


find_CCA_genes <- function(SeuratObjectList, num_hvg_genes=1000, num_samples_variable=2) {
  ob.list <- SeuratObjectList
  CCA_genes <- c()
  for (i in 1:length(ob.list)) {
    CCA_genes <- c(CCA_genes, head(rownames(ob.list[[i]]@hvg.info), num_hvg_genes))
  }
  CCA_genes <- names(which(table(CCA_genes) >= num_samples_variable))
  for (i in 1:length(ob.list)) {
    CCA_genes <- CCA_genes[CCA_genes %in% rownames(ob.list[[i]]@raw.data)]
  }
  return(CCA_genes)
}