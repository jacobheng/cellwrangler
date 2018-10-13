#'Create a single Soup Channel 
#'
#' createSoupChannel is a modified version of the SoupChannel() function in the excellent SoupX package. Instead 
#' of using a table of droplets (tod) that include all droplets to estimate the "soup", createSoupChannel feeds 
#' an expression matrix consisting of only empty droplets, determined by the testDrops() function, to estimate
#' the soup in the sample.
#'
#' @param blankDrops an expression matrix of empty droplets derived from the original raw matrix using 
#' the testDrops() function. Expression values from the blankDrops matrix are used to estimate the soup using
#' the emptySoup() function.
#' @param cells an expression matrix of cells (e.g. determined by the testDrops() function)
#' @param channelName string; the name of the channel 
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @param keepDroplets logical; whether to keep the blankDrops matrix after estimating the soup; defaults to
#' FALSE
#' @keywords createSoupChannel
#' @export
#' @return a SoupChannel object 
#' @examples
#' scl <- createSoupChannel(blankDrops=raw_mtx[,blank_barcodes], cells=raw_mtx[,cell_barcodes], channelNames=NULL
#' soupRange=c(0,10))

createSoupChannel <- function (blankDrops, cells, channelName, soupRange = c(0, 10), keepDroplets = FALSE, 
                               ...) 
{
  if (missing(channelName)) 
    channelName = "UnknownChannel"
  if (!all(gsub("___.*", "", colnames(blankDrops)) == channelName)) 
    colnames(blankDrops) = paste0(channelName, "___", colnames(blankDrops))
  if (!all(gsub("___.*", "", colnames(cells)) == channelName)) 
    colnames(cells) = paste0(channelName, "___", colnames(cells))
  out = list(blankDrops = blankDrops, toc = cells, channelName = channelName)
  out = c(out, list(...))
  out$nUMIs = Matrix::colSums(cells)
  out$nDropUMIs = Matrix::colSums(blankDrops)
  class(out) = c("list", "SoupChannel")
  channel = emptySoup(out, soupRange = soupRange, keepDroplets = keepDroplets)
  return(channel)
}