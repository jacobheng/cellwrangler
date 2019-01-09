#'Get genes for a kegg pathway
#'
#' get_kegg_pathway_genes() is a wrapper function around keggGet() that retrieves gene names in a KEGG pathway.
#'
#' @param pathway_name a KEGG pathway name
#' @param kegg_pathway_df a kegg pathway database that is produced using the cellwrangler get_kegg_database()
#' function using "pathway" as an option.
#' @keywords get_kegg_database
#' @export
#' @return a list
#' @examples
#' pathway_genes <- get_kegg_pathway_genes("Pentose phosphate pathway", mouse_kegg_pathways)

get_kegg_pathway_genes <- function(pathway_name, kegg_pathway_df){
  
  kegg_path_genes <- KEGGREST::keggGet(kegg_pathway_df[kegg_pathway_df$pathway %in% pathway_name,]$id)[[1]]$GENE
  kegg_path_genes <- cbind(kegg_path_genes)
  kegg_path_genes <- stringr::str_split_fixed(kegg_path_genes[grep(";", kegg_path_genes),], ";", n=2)[,1]
  kegg_path_genes <- as.vector(kegg_path_genes)
  return(kegg_path_genes)
  
}