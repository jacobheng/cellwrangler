#' Fit model function in monocle
#'
#' @description monocle_fitModel() is the fitModel() function found in monocle 2. This function fits a vector 
#' generalized additive model (VGAM) from the VGAM package for each gene in a CellDataSet. By default, expression 
#' levels are modeled as smooth functions of the Pseudotime value of each cell. That is, expression is a function 
#' of progress through the biological process. More complicated formulae can be provided to account for additional 
#' covariates (e.g. day collected, genotype of cells, media conditions, etc).
#' @param cds a monocle CellDataSet object
#' @param modelFormulaStr model formula string specifying the model to fit for the genes
#' @param relative_expr Whether to fit a model to relative or absolute expression. Only meaningful for count-based 
#' expression data. If TRUE, counts are normalized by Size_Factor prior to fitting.
#' @param cores the number of processor cores to be used during fitting.
#' @keywords monocle_fitModel() 
#' @export
#' @return a vglm object
#' @examples
#' monocle_fitModel(dat, modelFormulaStr="~num_genes_expressed+genotype")



monocle_fitModel <- function(cds, modelFormulaStr = "~sm.ns(Pseudotime, df=3)", 
          relative_expr = TRUE, cores = 1) 
{
  if (cores > 1) {
    f <- monocle:::mcesApply(cds, 1, monocle:::fit_model_helper, required_packages = c("BiocGenerics", 
                                                                   "Biobase", "VGAM", "plyr", "Matrix"), cores = cores, 
                   modelFormulaStr = modelFormulaStr, expressionFamily = cds@expressionFamily, 
                   relative_expr = relative_expr, disp_func = cds@dispFitInfo[["blind"]]$disp_func)
    f
  }
  else {
    f <- monocle:::smartEsApply(cds, 1, monocle:::fit_model_helper, convert_to_dense = TRUE, 
                      modelFormulaStr = modelFormulaStr, expressionFamily = cds@expressionFamily, 
                      relative_expr = relative_expr, disp_func = cds@dispFitInfo[["blind"]]$disp_func)
    f
  }
}