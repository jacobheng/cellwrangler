#'Get gene-level summary statistics for CellDataSet object
#'
#' summarize_cds_exprs() takes a monocle CellDataSet object as input and obtain summary statistics for each gene
#' for specified groups in the dataset. 
#'
#' @param cds a CellDataSet object 
#' @param group one or more column(s) specifying groups to obtain summary statistics for e.g. CellType
#' @param gene_short_name logical; whether to attach gene_short_name as a column to output dataframe
#' 
#' @keywords summarize_cds_exprs()
#' @export
#' @return a dataframe
#' @examples
#' dat_summary <- summarize_cds_exprs(dat, c("CellType"))

summarize_cds_exprs <- function(cds, group, gene_short_name = T) {
  
  if(length(group) == 1) {
    pData(cds)$group_variable <- pData(cds)[,group]
    
  } else {
    pData(cds)$group_variable <- interaction(pData(cds)[,group], sep= " - ")
  }
  group_vector <- as.character(unique(pData(cds)$group_variable))
  
  #Calculate means
  group_means <- sapply(group_vector, function(x) {
    means <- Matrix::rowMeans(exprs(cds[, pData(cds)$group_variable %in% x]))
    return(means)
  })
  
  means_melt <- melt(group_means)
  colnames(means_melt) <- c("gene_id", "group_variable", "mean")
  rownames(means_melt) <- interaction(means_melt$gene_id, means_melt$group_variable)
  
  #Calculate SD
  group_SD <- sapply(group_vector, function(x) {
    SD <- rowSD(exprs(cds[, pData(cds)$group_variable %in% x]))
    return(SD)
  })
  
  SD_melt <- melt(group_SD)
  colnames(SD_melt) <- c("gene_id", "group_variable", "std_dev")
  rownames(SD_melt) <- interaction(SD_melt$gene_id, SD_melt$group_variable)
  
  #Calculate SE
  group_SE <- sapply(group_vector, function(x) {
    SE <- rowSE(exprs(cds[, pData(cds)$group_variable %in% x]))
    return(SE)
  })
  SE_melt <- melt(group_SE)
  colnames(SE_melt) <- c("gene_id", "group_variable", "std_error")
  rownames(SE_melt) <- interaction(SE_melt$gene_id, SE_melt$group_variable)
  #Merge all
  stats_summary <- merge_df_list(list(means_melt, SD_melt, SE_melt), all.x=T, all.y=T)
  stats_summary <- stats_summary[,c("gene_id.x", "group_variable.x", "mean", "std_dev", "std_error")]
  if(length(group) == 1) {
    colnames(stats_summary) <- c("gene_id", group, "mean", "std_dev", "std_error")
  } else {
    colnames(stats_summary) <- c("gene_id", "group_variable", "mean", "std_dev", "std_error")
    n_groups <- length(group)
    for(i in 1:n_groups) {
      stats_summary[, group[[i]]] <- stringr::str_split_fixed(stats_summary$group_variable, pattern = " - ", 
                                                              n = n_groups )[,i]
    }
  }
  if(gene_short_name == T) {
  stats_summary$gene_short_name <- findGeneName(stats_summary$gene_id, cds, unique = F) 
  } else { stats_summary <- stats_summary}
  
  return(stats_summary)
}
