#' Find genes with high biological coefficent of variation
#'
#' find_highbcv_genes() return genes with a high biological coefficent of variation (bcv) from a monocle 
#' CellDataSet object, subject to a number of thresholds including expression in number of cells, expression
#' levels and residuals.
#' @param cds a monocle CellDataSet object
#' @param annotated_fData annotated fData object produced by the cellwrangler annotate_bcv() function; if NULL,
#' parameters of bcv will be calcuated from cds provided.
#' @param cell_threshold threshold for expression in number of cells
#' @param exprs_threshold threshold for expression levels
#' @param resid_threshold threshold for residuals
#' @keywords find_highbcv_genes
#' @export
#' @examples
#' find_highbcv_genes(dat)


find_highbcv_genes <- function(cds, annotated_fData = NULL, cell_threshold=10, exprs_threshold=0.001, 
                               resid_threshold=0.2) {
  if(is.null(annotated_fData) == T) {
  cds <- cds
  col_sums<-Matrix::colSums(exprs(cds))
  scaled_mtx <- t(t(exprs(cds))/col_sums)*10000
  fData(cds)$mean_exprs < -Matrix::rowMeans(scaled_mtx)
  cds_subsets<- lapply(seq(1,ceiling(nrow(dat0)/10000),1), function(x){
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
  fData(cds)$num_cells_expressed<-Matrix::rowSums(scaled_mtx>0.01)
  fData(cds)$percent_detection<-(fData(cds)$num_cells_expressed/dim(cds)[2])*100
  annotated_fData <- fData(cds)
  } else { annotated_fData <- annotated_fData }
  #fit model
  cds_fit<-mgcv::gam(log(bcv)~s(log(mean_exprs)),data=annotated_fData)
  cds_expressed_genes<-rownames(subset(annotated_fData,num_cells_expressed >= cell_threshold & mean_exprs >= exprs_threshold))
  high_bcv_genes_idx<-resid(cds_fit)>=resid_threshold
  high_bcv_genes<-intersect(rownames(resid(cds_fit))[high_bcv_genes_idx],cds_expressed_genes)
  return(high_bcv_genes)
}