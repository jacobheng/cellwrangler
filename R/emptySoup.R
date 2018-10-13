#'Estimate expression profile of soup from empty droplets
#'
#' emptySoup() is a modified version of the estimateSoup() function in the excellent SoupX package. Instead
#' of using all droplets within the specified soupRange to estimate the expression profile of the
#' soup as per the estimateSoup() function, emptySoup() uses only empty droplets e.g. determined using 
#' the testDrops() function to estimate the expression profile of the soup.
#'
#' @param sc a SoupChannel object created using the createSoupChannel() function implemented in cellwrangler.
#' @param soupRange Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.
#' @param keepDroplets logical; whether to keep the blankDrops matrix after estimating the soup; defaults to
#' FALSE
#' @keywords emptySoup
#' @export
#' @return a SoupChannel object 
#' @examples
#' scl <- emptySoup(scl, soupRange = c(0, 10), keepDroplets = FALSE)

emptySoup <- function (sc, soupRange = c(0, 10), keepDroplets = FALSE) 
{
  if (!is(sc, "SoupChannel")) 
    stop("sc must be a SoupChannel object.")
  w = which(sc$nDropUMIs > soupRange[1] & sc$nDropUMIs < soupRange[2])
  sc$soupProfile = SoupX:::estRateLims(Matrix::rowSums(sc$blankDrops[, w, drop = FALSE]), 
                                       sum(sc$blankDrops[, w]))
  if (!keepDroplets) 
    sc$blankDrops = NULL
  return(sc)
}