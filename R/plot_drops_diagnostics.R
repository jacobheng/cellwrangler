#'Plot diagnostic plot for testDrops results 
#'
#' @description Plots a diagnostic plot for a testDrops result. The x-axis represents the total number of UMIs 
#' per cell, while the y-axis represents the negative log probability that a droplet is a cell.
#'
#' @param testDrops_res a testDrops object create by the testDrops function
#' @param testDrops_res_list a list of testDrops objects create by the testDrops function; if a testDrops_res_list 
#' is supplied,testDrops_res should be the name of the testDrops object in that list i.e. testDrops_res = testDrops_res_list[[testDrops_res]]
#' @param x_range a vector specifying range of values to plot for x-axis
#' @param y_range a vector specifying range of values to plot for y-axis
#' @param log_x logical; whether x-axis should be log-transformed
#' @param log_y logical; whether y-axis should be log-transformed
#' @keywords plot_drops_diagnostics
#' @export
#' @return a ggplot2 object
#' @examples
#' plot_drops_diagnostics(testDrops_res)

plot_drops_diagnostics <- function(testDrops_res, testDrops_res_list=NULL, x_range=NULL, y_range=NULL, log_x=F, log_y=F  ) {
  if(is.null(testDrops_res_list) == F) {
    tmp <- as.data.frame(testDrops_res_list[[testDrops_res]]$data)
  } else { 
    tmp <- as.data.frame(testDrops_res$data)
    }
  
  p <- ggplot(data=tmp, aes(x=Total, y=-`LogProb`)) + geom_point(color=ifelse(tmp$is_cell == T, "red", "black")) + xlab("Total UMI count") + ylab("-Log Probability") + ggtitle(sample)
  p <- p + theme_classic() 
  
  if(is.null(x_range) == F) {
    p <- p + xlim(x_range[1], x_range[2])
  } else { p <- p }
  
  if(is.null(y_range) == F) {
    p <- p + ylim(y_range[1], y_range[2])
  } else { p <- p }
  
  if(log_x == T) {
    p <- p + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x)))
  } else { p <- p }
  
  if(log_y == T) {
    p <- p + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x)))
  } else { p <- p }
    
  return(p)
}


