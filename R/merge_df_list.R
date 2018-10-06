#' Merge a list of dataframes or matrices
#'
#' This is a wrapper function based on the merge_by_rownames function in cellwrangler for merging multiple
#' dataframes or matries by their rownames
#' @param df_list a list of dataframes or matrices
#' @param all.x logical;  if TRUE, then extra rows will be added to the output, one for each row in x that has no 
#' matching row in y. These rows will have NAs in those columns that are usually filled with values from y. 
#' Defaults to TRUE
#' @param all.y logical; analogous to all.x
#' @param sort logical; if TRUE, result will be sorted by columns. Defaults to FALSE
#' @keywords merge_df_list
#' @export
#' @return A merged dataframe or matrix with the appropriate rownames
#' @examples
#' dataframes_list <- list(dataframe1,dataframe2,dataframe3,dataframe4)
#' merge_df_list(dataframes_list)

merge_df_list <- function(df_list, all.x= TRUE, all.y = TRUE, sort = FALSE) {
  merged_df <- df_list[[1]]
  for(i in 2:length(df_list)) {
    merged_df <- merge_by_rownames(merged_df, df_list[[i]], all.x = all.x, all.y = all.y, sort = sort)
  }
  return(merged_df)
}