#'Plot biplot for sparse pca
#'
#' sparse_pca_biplot() creates a biplot for two principal components computed by sparse PCA (implemented 
#' by cellrangerRkit. 
#' 
#' @param sparse_pca_obj a sparse pca object created by the cellRangerRkit sparse_pca() function
#' @param pcs a two-component vector of PCs to plot e.g. c(1,2)
#' @param group a vector specifying groups for the observations. The points on the biplot will be colored
#' according to group
#' @keywords sparse_pca_biplot()
#' @export
#' @return A ggplot object
#' @examples
#' sparse_pca_biplot(sparse_pca_obj, pcs=c(1,2), group=pData(dat)$CellType)

sparse_pca_biplot <- function(sparse_pca_obj, pcs, group = NULL) {
  
  tmp <- as.data.frame(sparse_pca_obj$x[,pcs])
  colnames(tmp) <- c(paste("PC",pcs[1],sep=""), paste("PC",pcs[2],sep=""))

  pca_var <- as.data.frame((sparse_pca_obj$var_pcs)*100)
  pca_var$PC <- as.numeric(rownames(pca_var))
  colnames(pca_var) <- c("Percent_var_explained", "PC")
  
  if(is.null(group) == F) {
    tmp$group <- group
    p <- ggplot(tmp, aes_string(x=paste("PC",pcs[1],sep=""), y=paste("PC",pcs[2],sep="")))
    p <- p + geom_point(aes(color=group)) 
  } else {
    p <- ggplot(tmp, aes_string(x=paste("PC",pcs[1],sep=""), y=paste("PC",pcs[2],sep="")))
    p <- p + geom_point() 
  }
    p <- p + theme_classic()  
    p <- p + xlab(paste("PC",pcs[1]," (", pca_var[pca_var$PC == pcs[1],]$Percent_var_explained, 
                        " % of variance explained)",sep=""))
    p <- p + ylab(paste("PC",pcs[2]," (", pca_var[pca_var$PC == pcs[2],]$Percent_var_explained, 
                        " % of variance explained)",sep="")) 
    
  return(p)
}