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

plot_drops_diagnostics <- function(testDrops_res) {
  tmp <- as.data.frame(testDrops_res$data)
  p <- ggplot2::ggplot(data=tmp, aes(x=Total, y=-LogProb)) + ggplot2::geom_point(color=ifelse(tmp$is.cell == T, "red", "black")) 
  + ggplot2::xlab("Total UMI count") + ggplot2::ylab("-Log Probability") + ggplot2::ggtitle(sample)
  p <- p + ggplot2::theme_classic() + ggplot2::scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                                           labels = scales::trans_format("log10", scales::math_format(10^.x)))
}


