#'Plot variance explained by sparse PCA
#'
#' plot_sparse_pca_var() plots the variance explained by sparse PCA (implemented by cellrangerRkit) as a
#' line graph. 
#' 
#' @param sparse_pca_obj a sparse pca object created by the cellRangerRkit sparse_pca() function
#' @param n_pcs number of PCs to plot; defaults to NULL; if null, all the PCs will be plotted
#' @keywords plot_sparse_pca_var
#' @export
#' @return A ggplot object
#' @examples
#' plot_sparse_pca_var(sparse_pca_obj, n_pcs = 40)


plot_sparse_pca_var <- function(sparse_pca_obj, n_pcs) {
  
  tmp <- as.data.frame((sparse_pca_obj$var_pcs)*100)
  tmp$PC <- as.numeric(rownames(tmp))
  colnames(tmp) <- c("Percent_var_explained", "PC")
  if(is.null(n_pcs) == F) {
    tmp <- tmp[n_pcs,]
  } else { tmp <- tmp }
  
  p <- ggplot(tmp, aes(x=PC, y=Percent_var_explained)) + geom_line() + geom_point() 
  p <= p + scale_x_continuous(breaks=tmp$PC) + xlab("PC") + ylab("Percentage of variance explained")
  
  return(p)
  
}