#'Load an expression matrix from a directory of cellranger-aligned 10X Genomics scRNAseq data 
#'
#' @description load_10X_matrix() loads an expression matrix as a sparseMatrix from a directory of 10X Genomics 
#' scRNAseq data created using the cellranger pipeline. Either the raw or the filtered matrix can be loaded.
#'
#' @param cellranger_outs_path the path to the "outs" directory in the cellranger library folder e.g.
#' "/mycellranger_library/outs"
#' @param cellranger_filter whether to use the filtered matrices from the cellranger pipeline as the 
#' expression matrix (i.e. as cells)
#' @param which_matrix must be "raw" or "filtered" to specify the raw or filtered matrix to load respectively.
#' @param cellranger_v3 logical; if cellranger version 3 and above was used to produced matrices.
#' @keywords load_10X_matrix
#' @export
#' @return a dgTMatrix 
#' @examples
#' raw_mtx <- load_10X_matrix("/mycellranger_library/outs", which_matrix="raw")


load_10X_matrix <- function(cellranger_outs_path, which_matrix="raw", cellranger_v3 = T) {
  
  file_names <- list.files(cellranger_outs_path, recursive = T)
  if(which_matrix == "filtered") {
    
    matrix_path <- file_names[grep("filtered.*matrix.mtx" ,file_names)]
    barcodes_path <- file_names[grep("filtered.*barcodes.tsv" ,file_names)]
    
    if(cellranger_v3 == T) { genes_path <- file_names[grep("filtered.*features.tsv" ,file_names)] 
    } else { genes_path <- file_names[grep("filtered.*genes.tsv" ,file_names)] }
    
  } else { 
    if(which_matrix == "raw") {
      
      matrix_path <- file_names[grep("raw.*matrix.mtx" ,file_names)]
      barcodes_path <- file_names[grep("raw.*barcodes.tsv" ,file_names)]
      
      if(cellranger_v3 == T) { genes_path <- file_names[grep("raw.*features.tsv" ,file_names)] 
      } else { genes_path <- file_names[grep("raw.*genes.tsv" ,file_names)] 
      }
      
    } else { message("Need to specify 'raw' or 'filtered' in which_matrix parameter") }
  }
  
  #Read in matrix
  mtx <- Matrix::readMM(paste0(cellranger_outs_path, matrix_path))
  
  #Append barcodes
  mtx_barcodes <- read.delim(paste0(cellranger_outs_path, barcodes_path), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
  colnames(mtx_barcodes) <- c("barcode")
  rownames(mtx_barcodes) <- mtx_barcodes$barcode
  colnames(mtx) <- mtx_barcodes$barcode
  
  #Read in genes/features
  mtx_genes <- read.delim(paste0(cellranger_outs_path, genes_path), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
  #Add colnames for fData
  if(cellranger_v3 == T) {  colnames(mtx_genes) <- c("id","gene_short_name", "feature_type")
  } else { 
    colnames(mtx_genes) <- c("id","gene_short_name") 
  }
  rownames(mtx_genes) <- mtx_genes$id
  rownames(mtx) <- mtx_genes$id
  
  return(mtx)
}