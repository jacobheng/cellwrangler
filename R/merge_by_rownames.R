
#' Merge dataframes or matrices by rownames
#'
#' This is a wrapper function based on the base R merge function for merging 
#' two dataframes or matrices by their rownames and returing a dataframe or matrix with the appropriate rownames.
#' Unlike merge, this function does not create an additional column for rownames.
#' @param x,y dataframes or matrices
#' @param all.x logical;  if TRUE, then extra rows will be added to the output, one for each row in x that has no 
#' matching row in y. These rows will have NAs in those columns that are usually filled with values from y. 
#' Defaults to TRUE
#' @param all.y logical; analogous to all.x
#' @param sort logical; if TRUE, result will be sorted by columns. Defaults to FALSE
#' @keywords merge.by.rownames
#' @export
#' @examples
#' merge_by_rownames(data.frame1, data.frame2)

merge_by_rownames <- function(x, y, all.x= TRUE, all.y = TRUE, sort = FALSE) {
  tmp <- merge(x, y, by = 0, all.x = all.x, all.y = all.y, sort = sort)
  rownames(tmp) <- tmp[,1]
  tmp <- tmp[,-1]
  return(tmp)
}