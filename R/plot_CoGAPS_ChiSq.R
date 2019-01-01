#'Plots the mean ChiSq values of a set of CoGAPS results
#'
#' plot_CoGAPS_ChiSq() plots the mean ChiSq values of a set of CoGAPS results.
#' 
#' @param CoGAPS_res_set a CoGAPS result set.
#' @keywords plot_CoGAPS_ChiSq
#' @export
#' @return A ggplot2 object
#' @examples
#' plot_CoGAPS_ChiSq(myCoGAPSres)


plot_CoGAPS_ChiSq <- function(CoGAPS_set) {
  nRuns <- length(CoGAPS_set)
  ChiSq <- sapply(seq(1,nRuns,1), function(x) {
    CoGAPS_set[[x]]@metadata$meanChiSq
  })
  tmp <- as.data.frame(names(CoGAPS_set))
  colnames(tmp) <- c("run")
  tmp$ChiSq <- ChiSq
  print(tmp)
  p <- ggplot(tmp, aes(x=run, y=ChiSq, group=1)) + geom_point() + geom_line()
  return(p)
}