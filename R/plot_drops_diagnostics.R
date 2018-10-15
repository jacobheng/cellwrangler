#'Plot diagnostic plot for testDrops results 
#'
#' Plots a diagnostic plot for a testDrops result. The x-axis represents the total number of UMIs per cell,
#' while the y-axis represents the negative log probability that a droplet is a cell.
#'
#' @param testDrops_res a testDrops object create by a testDrops function
#' @keywords plot_drops_diagnostics
#' @export
#' @return a ggplot2 object
#' @examples
#' plot_drops_diagnostics(testDrops_res)

plot_drops_diagnostics <- function(sample, testDrops_res_list) {
  tmp <- as.data.frame(testDrops_res_list[[sample]]$data)
  p <- ggplot(data=tmp, aes(x=Total, y=-LogProb)) + geom_point(color=ifelse(tmp$is.cell == T, "red", "black")) 
  + xlab("Total UMI count") + ylab("-Log Probability") + ggtitle(sample)
  p <- p + theme_classic() + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                           labels = scales::trans_format("log10", scales::math_format(10^.x)))
}


