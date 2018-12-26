#' Fit vector generalized linear model to a CellDataSet object
#'
#' fit_cds_vglm() is a wrapper function around the monocle::fitModel() function that fits a vglm for each 
#' gene in a CellDataSet object and return a dataframe with the estimate, standard error, z value and
#' p-value for each coefficient specified in the model formula string.
#' @param cds a monocle CellDataSet object
#' @param modelFormulaStr model formula string specifying the model to fit for the genes
#' @param test_genes genes to test the model on; test_genes must be expressed in cds object.
#' @keywords fit_cds_vglm()
#' @export
#' @examples
#' find_cds_vglm(dat)

fit.cds.vglm <- function(cds, modelFormulaStr, test_genes=c("Actb")){
  #Test function
  test.VGAM <- monocle::fitModel(cds[findGeneID(test_genes, cds),], modelFormulaStr = modelFormulaStr)
  test.VGAM.summary <- summaryvglm(test.VGAM[[1]])
  #Error trapping
  robust_summaryvglm <- function(x) {tryCatch(summaryvglm(x),warning = function(w) {print(paste(w))},  
                                              error = function(e) {print(paste(e))}) }
  
  robust_coef <- function(x) {tryCatch(coef(x), error = function(e) { 
    tmp <- matrix(0,nrow=nrow(coef(test.VGAM.summary)),ncol=ncol(coef(test.VGAM.summary)))
    colnames(tmp) <- colnames(coef(test.VGAM.summary))
    rownames(tmp) <- rownames(coef(test.VGAM.summary))
    return(tmp)
  }) }
  #Fit vglm
  cds.vglm.list <- monocle::fitModel(cds, modelFormulaStr = modelFormulaStr)
  cds.vglm_coef <- lapply(names(cds.vglm.list),function(x){
    tmp1 <- robust_summaryvglm(cds.vglm.list[[x]])
    tmp2 <- melt(robust_coef(tmp1))
    tmp2$interaction <- interaction(tmp2[,1:2],sep=":")
    tmp2 <- t(tmp2[,c("interaction","value")])
    colnames(tmp2) <- tmp2[1,]
    tmp2 <- t(as.data.frame(tmp2[-1,]))
    rownames(tmp2) <- x
    return(tmp2)
  })
  ncol_coef <- (length(unlist(str_split(modelFormulaStr,pattern="\\+")))+1)*4
  cds.vglm_coef <- cds.vglm_coef[unlist(lapply(cds.vglm_coef, function(x){ 
    ncol(x) == ncol_coef
  }))]
  cds.vglm_coef.df <- as.data.frame(do.call(rbind, cds.vglm_coef))
  coef.cols <- colnames(cds.vglm_coef.df)
  cds.vglm_coef.df[coef.cols] <- lapply(cds.vglm_coef.df[coef.cols], function(x){
    num.x <- as.numeric(as.character(x))
    return(num.x)
  })
  cds.vglm.merged <- as.data.frame(merge.by.rownames(cds.vglm_coef.df,fData(cds)))
  return(cds.vglm.merged)
}