#'Load a set of CoGAPS results
#'
#' load_CoGAPS_res() loads a set of CoGAPS results.
#'
#' @param res_dir name of directory containing CoGAPS results. Each CoGAPS result should be contained in
#' a folder with the name "run_x" in numerical order baesed on the length of nFactor_range.
#' @param res_name name of .RData result files. 
#' @param nFactor_range range of nFactor used in the set of CoGAPS results.
#' 
#' @keywords load_CoGAPS_res
#' @export
#' @return a list of CoGAPS objects
#' @examples
#' myCoGAPSres <- load_CoGAPS_res(res_dir = "myCoGAPSres_set", res_name = "myCoGAPS_result", 
#' nFactor_range = c(25,60,5))

load_CoGAPS_res <- function(res_dir, res_name, nFactor_range) {
  
  nFactor_range <- nFactor_range
  nPatterns <- lapply(nFactor_range,function(x){
    paste("nP",x,sep="")
  })
  
  CoGAPS_res <- lapply(seq(1,length(nPatterns),1), function(x){
    get(load(paste(res_dir,"/run_",x, "/", res_name,".RData",sep=""), verbose=T))
  })
  names(CoGAPS_res) <- nPatterns
  return(CoGAPS_res)
}






