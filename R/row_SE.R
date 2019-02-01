#'Get standard error of the mean for rows in a matrix
#'
#' rowSE() calculates row-wise standard error of the mean for a matrix
#'
#' @param matrix a matrix
#' 
#' @keywords rowSE()
#' @export
#' @return a numerical vector
#' @examples
#' row_Std_Err <- rowSE(matrix)

rowSE <- function(matrix) {
  std_error <- function(x) { sd(x)/sqrt(length(x)) }
  SE <- apply(matrix, 1, std_error)
  return(SE)
}