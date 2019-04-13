#'Get standard deviation for rows in a matrix
#'
#' @description rowSD() calculates row-wise standard deviation for a matrix
#'
#' @param matrix a matrix
#' 
#' @keywords rowSD()
#' @export
#' @return a numerical vector
#' @examples
#' row_Std_Dev <- rowSD(matrix)

rowSD <- function(matrix) {
  SD <- apply(matrix, 1, sd)
  return(SD)
}
