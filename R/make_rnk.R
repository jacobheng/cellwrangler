#'Make rank lists
#'
#' @description make_rnk() makes rank lists for pre-ranked geneset enrichment analysis.
#'
#' @param celltype CellType to create rank lists for.
#' @param rank_df a dataframe with the pre-requisite rank statistics
#' @param gene_name column in rank_df containing gene names.
#' @param rank_stat column in rank_df containing statistic to rank genes by.
#' @param group_column column name of column in annotation object with values corresponding to those 
#' in group_vector to use for computing logical columns. Defaults to NULL
#' @param decreasing logical; whether to rank genes by decreasing value of rank_stat.
#' @param write_rnk logical; whether to write .rnk files for GSEA
#' @param output_dir name of directory to write .rnk files to; if directory does not exist, it will be created
#' @param output_name prefix for .rnk file names; .rnk files will be named "output_name_celltype.rnk"
#' @keywords make_rnk
#' @export
#' @return a list
#' @examples
#' make_rnk(Rods, rank_df, "Gene_symbol", "Rank stat", decreasing =T, write_rnk=F)


make_rnk <- function(celltype, rank_df, gene_name, rank_stat, decreasing = T, 
                     write_rnk=F, output_dir=NULL, output_name=NULL) {
  tmp <-  rank_df[rank_df$myCellType == celltype, c(gene_name,rank_stat)]
  print(paste("Ranking", celltype, gene_name, "by", rank_stat, sep = " "))
  tmp <- tmp[order(tmp[,rank_stat], decreasing = decreasing),]
  tmp <- tmp[!tmp[,gene_name] == "#N/A",]
  tmp <- tmp[!duplicated(tmp[,gene_name]),]
  tmp <- tmp[!is.na(tmp[,gene_name]),]
  rownames(tmp) <- as.character(tmp[,gene_name])
  tmp_vector <- tmp[,2]
  names(tmp_vector) <- rownames(tmp)
  
  if(write_rnk == T) {
    if(dir.exists(output_dir) == F ) {
      dir.create(output_dir)
    } else {    }
  write.table(tmp, paste(output_dir,"/", output_name, "_",celltype,".rnk", sep=""), row.names = F, 
              quote = F, sep = "\t", eol = "\r\n")
  print( paste(output_dir,"/", output_name, "_",celltype,".rnk saved", sep="") )
  } else { tmp <- tmp }
  
  return(tmp_vector)
}