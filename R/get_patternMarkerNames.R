#'Get pattern markers as gene names for a set of CoGAPS results
#'
#' @description get_patternMarkerNames is a wrapper function around the CoGAPS::patternMarkers function for 
#' retrieving pattern markers as gene names for a set of CoGAPS results derived from a monocle CellDataSet object.
#'
#' @param CoGAPS_res a CoGAPS result object.
#' @param cds a CellDataSet object e.g. used in the monocle package
#' @keywords get_patternMarkerNames
#' @export
#' @return a CoGAPS patternMarkers object
#' @examples
#' patternMarkers <- get_patternMarkerNames(myCoGAPSres, cds=dat)



get_patternMarkerNames <- function(CoGAPS_res, cds, threshold = "all", lp = NA) {
  tmp <- CoGAPS_res
  sd_rownames <- findGeneName(rownames(CoGAPS_res@featureStdDev), cds)
  rownames(tmp@featureStdDev) <- make.names(sd_rownames, unique = T)
  loading_rownames <- findGeneName(rownames(CoGAPS_res@featureLoadings))
  rownames(tmp@featureLoadings) <- make.names(loading_rownames, unique = T)
  pM <- patternMarkers(tmp, threshold = threshold, lp = lp)
  return(pM)
}