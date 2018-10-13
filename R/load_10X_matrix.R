#'Load an expression matrix from a directory of cellranger-aligned 10X Genomics scRNAseq data 
#'
#' load_10X_matrix() loads an expression matrix as a sparseMatrix from a directory of 10X Genomics scRNAseq 
#' data created using the cellranger pipeline. Either the raw or the filtered matrix can be loaded.
#'
#' @param cellranger_outs_path the path to the "outs" directory in the cellranger library folder e.g.
#' "/mycellranger_library/outs"
#' @param cellranger_filter whether to use the filtered matrices from the cellranger pipeline as the 
#' expression matrix (i.e. as cells)
#' @param which_matrix must be "raw" or "filtered" to specify the raw or filtered matrix to load respectively.
#' @keywords load_10X_matrix
#' @export
#' @return a dgTMatrix 
#' @examples
#' raw_mtx <- load_10X_matrix("/mycellranger_library/outs", which_matrix="raw")


load_10X_matrix <- function(cellranger_outs_path, which_matrix="raw") {
  genome <- list.files(file.path(cellranger_outs_path, "raw_gene_bc_matrices_mex"))
  if(which_matrix == "filtered"){
    #Read in filtered mtx
    mtx <- Matrix::readMM(paste(cellranger_outs_path,"/filtered_gene_bc_matrices_mex/", genome,"/matrix.mtx", sep=""))
    #Filtered mtx genes
    mtx_genes <- read.delim(paste(cellranger_outs_path,"/filtered_gene_bc_matrices_mex/", genome,"/genes.tsv", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
    colnames(mtx_genes) <- c("id","gene_short_name")
    rownames(mtx_genes) <- mtx_genes$id
    rownames(mtx) <- mtx_genes$id
    #Append barcodes
    mtx_barcodes <- read.delim(paste(cellranger_outs_path,"/filtered_gene_bc_matrices_mex/", genome,"/barcodes.tsv", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
    colnames(mtx_barcodes) <- c("barcode")
    rownames(mtx_barcodes) <- mtx_barcodes$barcode
    colnames(mtx) <- mtx_barcodes$barcode
  } else { if(which_matrix=="raw") {
    #Read in raw mtx
    mtx <- Matrix::readMM(paste(cellranger_outs_path,"/raw_gene_bc_matrices_mex/", genome,"/matrix.mtx", sep=""))
    #Raw mtx genes
    mtx_genes <- read.delim(paste(cellranger_outs_path,"/raw_gene_bc_matrices_mex/", genome,"/genes.tsv", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
    colnames(mtx_genes) <- c("id","gene_short_name")
    rownames(mtx_genes) <- mtx_genes$id
    rownames(mtx) <- mtx_genes$id
    #Append barcodes
    mtx_barcodes <- read.delim(paste(cellranger_outs_path,"/raw_gene_bc_matrices_mex/", genome,"/barcodes.tsv", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
    colnames(mtx_barcodes) <- c("barcode")
    rownames(mtx_barcodes) <- mtx_barcodes$barcode
    colnames(mtx) <- mtx_barcodes$barcode
  } else { message("Need to specify 'raw' or 'filtered' in which_matrix parameter") }
  } 
  return(mtx)
}