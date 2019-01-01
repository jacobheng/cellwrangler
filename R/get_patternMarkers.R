#'Get pattern markers for a set of CoGAPS results
#'
#' get_patternMarkers is a wrapper function around the CoGAPS::patternMarkers function for retrieving pattern 
#' markers for a set of CoGAPS results.
#'
#' @param CoGAPS_res_set a CoGAPS result set. If set to NULL, a Pmeans matrix may be supplied directly as input
#' for the Pattern_set parameter.
#' @param nFactor_range range of nFactor used in the set of CoGAPS results.
#' @param relabel_genes logical; if genes should be relabeled from Ensembl IDs to gene symbols.
#' @param cds CellDataSet object for relabeling genes.
#' 
#' @keywords load_CoGAPS_res
#' @export
#' @return a list of CoGAPS objects
#' @examples
#' myCoGAPSres <- load_CoGAPS_res(res_dir = "myCoGAPSres_set", res_name = "myCoGAPS_result", 
#' nFactor_range = c(25,60,5))



get_patternMarkers <- function(CoGAPS_res, nFactor_range, relabel_genes=F, cds=NULL,
                               scaledPmatrix=T, full=T) {
  nFactor_range <- nFactor_range
  nPatterns <- lapply(nFactor_range,function(x){
    paste("nP",x,sep="")
  })
    
    pM <- lapply(nPatterns, function(n){
      labeled_Amatrix <- CoGAPS_res[[n]]$Amean
      if(relabel_genes == T) {
        rownames(labeled_Amatrix) <- findGeneName(cds, rownames(labeled_Amatrix))
      } else { labeled_Amatrix <- labeled_Amatrix  }

      tmp <- patternMarkers(Amatrix=labeled_Amatrix, scaledPmatrix=scaledPmatrix, full=full) 
      #patternMarkers return list of vectors
      #make each column in tmp the same length, and rbind into dataframe
      tmp$PatternMarkers <- t(do.call(rbind,lapply(tmp$PatternMarkers, function(x) c(x,rep(NA,max(sapply(tmp$PatternMarkers,length)-length(x)))))))
      colnames(tmp$PatternMarkers) <- seq(1:ncol(tmp$PatternMarkers))
      return(tmp)
    })
    names(pM) <- nPatterns
  
  return(pM)
}