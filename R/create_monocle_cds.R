#'Create a monocle CellDataSet object from a directory of cellranger-aligned 10X Genomics scRNAseq data 
#'
#' create_monocle_cds() creates a monocle CellDataSet object from a directory of 10X Genomics scRNAseq data
#' created using the cellranger pipeline. Includes options to determine which barcodes are cells vs non-cells
#' in the original raw matrix.
#'
#' @param cellranger_outs_path the path to the "outs" directory in the cellranger library folder e.g.
#' "/mycellranger_library/outs"
#' @param which_matrix must be "raw" or "filtered" to specify the raw or filtered matrix to load respectively.
#' @param cellranger_v3 logical; if cellranger version 3 and above was used to produced matrices.
#' @param UMI_cutoff numeric; a single cutoff of total UMI counts per cell to use for determining cells (i.e. 
#' all barcodes with total UMI counts above this cutoff will be included in the expression matrix); defaults 
#' to NULL
#' @param testDrops logical; whether to use the testDrops function to determine which barcodes are cells (and
#' therefore included in the expression matrix)
#' @param lower_UMI_threshold numeric; parameter for the testDrops function. All barcodes with total UMI counts
#' below this value are considered empty droplets (i.e. not cells)
#' @param upper_UMI_threshold numeric; parameter for the testDrops function. All barcodes with total UMI counts
#' above this value are considered cells
#' @param FDR numeric; parameter for the testDrops function. FDR value to use in the testDrops function when
#' statistically determining cells vs non-cells
#' @param cell_barcodes a vector of barcodes to use for creating the expression matrix (all barcodes specified
#' by cell_barcodes are assumed to be cells)
#' @keywords create_monocle_cds
#' @export
#' @return a monocle CellDataSet object
#' @examples
#' dat <- create_monocle_cds(cellranger_outs_path, cellranger_filter= F, UMI_cutoff=NULL, testDrops=T,
#' lower_UMI_threshold=100, upper_UMI_threshold=500, FDR=0.01, cell_barcodes=NULL)

create_monocle_cds <- function(cellranger_outs_path, cellranger_v3 = T,  which_matrix="raw", UMI_cutoff=NULL, testDrops=F,
                               lower_UMI_threshold=100, upper_UMI_threshold=500, FDR=0.01, cell_barcodes=NULL) {
  
  if(cellranger_v3 == T) {
    if(which_matrix == "filtered"){
      #Read in filtered mtx
      mtx <- Matrix::readMM(paste(cellranger_outs_path,"/filtered_feature_bc_matrix/matrix.mtx.gz", sep=""))
      #Filtered mtx genes
      mtx_genes <- read.delim(paste(cellranger_outs_path, "/filtered_feature_bc_matrix/features.tsv.gz", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
      if(unique(mtx_genes[,3]) == "Gene Expression") {
        colnames(mtx_genes) <- c("id","gene_short_name", "feature_type")
        rownames(mtx_genes) <- mtx_genes$id
        rownames(mtx) <- mtx_genes$id } else { stop("Only gene expression data are supported by cellwrangler") }
      #Append barcodes
      mtx_barcodes <- read.delim(paste(cellranger_outs_path,"/filtered_feature_bc_matrix/barcodes.tsv.gz", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
      colnames(mtx_barcodes) <- c("barcode")
      rownames(mtx_barcodes) <- mtx_barcodes$barcode
      colnames(mtx) <- mtx_barcodes$barcode
    } else { if(which_matrix=="raw") {
      #Read in raw mtx
      mtx <- Matrix::readMM(paste(cellranger_outs_path,"/raw_feature_bc_matrix/matrix.mtx.gz", sep=""))
      #Raw mtx genes
      mtx_genes <- read.delim(paste(cellranger_outs_path,"/raw_feature_bc_matrix/features.tsv.gz", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
      if(unique(mtx_genes[,3]) == "Gene Expression") {
        colnames(mtx_genes) <- c("id","gene_short_name", "feature_type")
        rownames(mtx_genes) <- mtx_genes$id
        rownames(mtx) <- mtx_genes$id } else { stop("Only gene expression data are supported by cellwrangler") }
      #Append barcodes
      mtx_barcodes <- read.delim(paste(cellranger_outs_path,"/raw_feature_bc_matrix/barcodes.tsv.gz", sep=""), stringsAsFactors = FALSE, sep = "\t", header = FALSE)
      colnames(mtx_barcodes) <- c("barcode")
      rownames(mtx_barcodes) <- mtx_barcodes$barcode
      colnames(mtx) <- mtx_barcodes$barcode
    } else { message("Need to specify 'raw' or 'filtered' in which_matrix parameter") }
    } 
    
  } else {
    
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
  }

  
  
  #Create monocle cds
  cds <- monocle::newCellDataSet(cellData=as(mtx, "sparseMatrix"),
                                 phenoData=new("AnnotatedDataFrame",data=mtx_barcodes),
                                 featureData=new("AnnotatedDataFrame",data=mtx_genes),
                                 lowerDetectionLimit=1,
                                 expressionFamily=VGAM::negbinomial.size())
  
  #Apply single cutoff
  if(is.null(UMI_cutoff) == F) {
    #Compute UMI totals for each cell
    UMI_count <- as.data.frame(Matrix::colSums(mtx))
    colnames(UMI_count) <- c("UMI_count")
    rownames(UMI_count) <- mtx_barcodes[,1]
    cutoff_barcodes <- rownames(UMI_count[UMI_count$UMI_count >= UMI_cutoff, , drop=F])
    cds <- cds[,cutoff_barcodes]
  } else { cds <- cds }
  #Determine cells with testDrops function
  if(testDrops==T) {
    mtx_split<- split_aggr_mtx(mtx)
    testDrops_res <- lapply(mtx_split, testDrops, lower_UMI_threshold=lower_UMI_threshold,
                            upper_UMI_threshold=upper_UMI_threshold, test.ambient=F, FDR=FDR)
    cell_drops<- c(unlist(lapply(sample_drops, function(x) {x$cell_barcodes})))
    cds <- cds[,cell_drops]
  } else { cds <- cds}
  #Select cell barcodes
  if(is.null(cell_barcodes) == F) {
    cds <- cds[,cell_barcodes]
  } else {cds <- cds}
  
  return(cds)
}