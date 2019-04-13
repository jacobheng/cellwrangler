#' Find rows by row names
#'
#' @description findRows() return specified rows in a matrix or dataframe by specifying the corresponding rownames.
#' @param row_names a vector of rowname(s) of row(s) to return
#' @param matrix matrix or dataframe
#' @keywords findRows
#' @export
#' @examples
#' findRows("myrowname", mydataframe)


findRows <- function(row_names, matrix) {
  tmp <- matrix[rownames(matrix) %in% row_names,]
  return(tmp)
}