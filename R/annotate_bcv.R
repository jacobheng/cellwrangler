#' Calculate and annotate biological coefficent of variation
#'
#' @description annotate_bcv() calculates the biological coefficent of variation from a give cds and returns an annotated
#' fData dataframe.
#' @param cds a monocle CellDataSet object
#' @param scaled_matrix a scaled expression matrix derived from cds provided. Exprs(cds) should be divided by 
#' the colsums and multiplied by 10,000.
#' @keywords annotate_bcv
#' @export
#' @return an annotated fData dataframe
#' @examples
#' annotate_bcv(dat)


annotate_bcv <- function(cds, scaled_matrix = NULL) {
  cds <- cds
  
  if(is.null(scaled_matrix) == F) {
    scaled_mtx <- scaled_matrix
  } else {
    col_sums<-Matrix::colSums(exprs(cds))
    scaled_mtx <- Matrix::t( Matrix::t( exprs(cds) )/col_sums)*10000
  }
  
  fData(cds)$mean_exprs<-Matrix::rowMeans(scaled_mtx)
  cds_subsets<- lapply(seq(1,ceiling(nrow(cds)/10000),1), function(x){
    tmp <- scaled_mtx[((x-1)*10000+1):min(c((x)*10000, max(nrow(cds)))),] 
    return(tmp) })
  subset_sd <- lapply(cds_subsets, function(x){
    tmp <- as.matrix(apply(x,1,sd))
    return(tmp)
  } )
  cds_sd <- do.call(rbind, subset_sd)
  print(table(rownames(cds_sd) == rownames(exprs(cds))))
  fData(cds)$sd_expr <- cds_sd
  fData(cds)$bcv<-(fData(cds)$sd_expr/fData(cds)$mean_exprs)**2
  fData(cds)$num_cells_expressed<-Matrix::rowSums(scaled_mtx> 0 )
  fData(cds)$percent_detection<-(fData(cds)$num_cells_expressed/dim(cds)[2])*100
  annotated_fData <- fData(cds)
  return(annotated_fData)
}