#'Get standard deviation or standard error for rows in a matrix
#'
#' rowSD() and rowSE() calculates row-wise standard deviation and standard error of the mean for a matrix
#'
#' @param matrix a matrix
#' 
#' @keywords rowSD(), rowSE()
#' @export
#' @return a numerical vector
#' @examples
#' row_Std_Dev <- rowSD(matrix)

rowSD <- function(matrix) {
  SD <- apply(matrix, 1, sd)
  return(SD)
}

rowSE <- function(matrix) {
  std_error <- function(x) { sd(x)/sqrt(length(x)) }
  SE <- apply(matrix, 1, std_error)
  return(SE)
}