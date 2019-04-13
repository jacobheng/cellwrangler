#'Get kegg database for an organism
#'
#' @description get_kegg_database() is a wrapper function around the keggList() function that retrieves 
#' a kegg database for an organism and returns a dataframe
#'
#' @param database a KEGG database; one of options available via KEGGREST listDatabases()
#' @param organism KEGG organism code e.g. "mmu" for "Mus musculus"
#' @keywords get_kegg_database
#' @export
#' @return a dataframe
#' @examples
#' kegg_db <- get_kegg_database("pathway","mmu")


get_kegg_database <- function(database, organism) {
  
  organism_ref <- as.data.frame((keggList("organism")))
  species <- organism_ref[organism_ref$organism == organism,]$species
  kegg_database <- KEGGREST::keggList(database, organism)
  kegg_database <- as.data.frame(cbind(kegg_database))
  kegg_database$id <- stringr::str_split_fixed(rownames(kegg_database), ":", n =2)[,2]
  kegg_database$species <- species
  species_tag <- stringr::str_split_fixed(species, " ", n=2)[,1]
  kegg_database[database] <- stringr::str_split_fixed(kegg_database[,1], paste(" - ", species_tag, sep=""), n = 2)[,1]
  
  return(kegg_database)
  
}