
#'Create multiple Soup Channels from an cellranger-aggregated matrix
#'
#' @description createMultiSoupChannel() is a wrapper function around the createSoupChannel() function that creates
#' a channel for each sample in a matrix aggregated from multiple samples with the cellranger pipeline.
#' The function requires a AllBlankDrops object, which is a matrix consisting of the empty droplets from
#' all samples derived from the original cellranger-aggregated raw matrix using the testDrops() function.
#' Based on the excellent SoupX package.
#' @param cds a CellDataSet object (e.g. used in the monocle package)
#' @param AllBlankDrops an expression matrix of empty droplets from all samples in the cellranger aggregate.
#' Matrix must contain barcodes from all the different samples used to form the cellranger aggregate; this
#' can be derived from the original aggregated raw matrix using the testDrops function. Expression values
#' from the AllBlankDrops matrix are used to estimate the soup.
#' @param channelNames vector of strings; channel names to use for all the channels in the sample; number
#' of channel names must equal number of samples used in the cellranger aggregate
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @keywords createMultiSoupChannel
#' @export
#' @return a SoupChannel object with multiple channels (corresponding to number of samples)
#' @examples
#' scl <- createMultiSoupChannel(cds=dat, AllBlankDrops=raw_mtx[,all_blank_barcodes], channelNames=NULL
#' soupRange=c(0,10))

createMultiSoupChannel <- function(cds, AllBlankDrops, channelNames = NULL, soupRange=c(0,10)) {
  if (is.null(channelNames)) {
    cell_barcodes <- rownames(pData(cds))
    channelNames <- sprintf("Channel%s", unique(str_split_fixed(cell_barcodes,pattern="-",n=2)[,2]))
  } else {
    channelNames = channelNames
  }
  channels = list()
  for (i in seq_along(channelNames)) {
    message(paste("Loading data for 10x",channelNames[i]))
    
    cell_matrix = exprs(cds)[,stringr::str_split_fixed(colnames(exprs(cds)),pattern="-",n=2)[,2]== i]
    blank_drops = AllBlankDrops[,stringr::str_split_fixed(colnames(AllBlankDrops),pattern="-",n=2)[,2]== i]
    
    channels[[channelNames[i]]] = createSoupChannel(blankDrops = blank_drops, cells = cell_matrix, 
                                                    channelName = channelNames[i], soupRange = soupRange, keepDroplets = F)
  }
  channels = SoupChannelList(channels)
  return(channels)
}
