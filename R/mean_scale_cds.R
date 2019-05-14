#'Mean scale expression of cell groups across multiple samples to have the same mean in CellDataSet object
#'
#' @description mean_scale_cds() scales the expression of cell groups e.g. CellType across multiple samples in
#' a monocle CellDataSet object such that cells from the same cell group in every sample have the same mean
#' expression. Assuming cells in same cell group have similar variances aross all samples, this will allow cells
#' from the same cell group in every sample to have similar distributions. 
#' @param cds monocle CellDataSet object
#' @param sample_col name of column in pData(cds) that annotates samples that cells are from
#' @param group name of column in pData(cds) that annotates groups that cells are in. 
#' @keywords mean_scale_cds
#' @export
#' @return A CellDataSet object
#' @examples
#' dat <- mean_scale_cds(dat0, sample_col="sample", group="CellType")

mean_scale_cds <- function (cds, sample_col, group) 
{
  cds <- cds
  group_values <- unique(pData(cds)[, group])
  samples <- unique(pData(cds)[, sample_col])
  print(samples)
  get_group_mean <- function(x, cds, y = NULL) {
    if (is.null(y) == TRUE) {
      tmp <- mean(Matrix::colSums(exprs(cds[, pData(cds)[,group] == 
                                              x])))
    }
    else {
      tmp <- mean(Matrix::colSums(exprs(cds[, pData(cds)[,group] == 
                                              x & pData(cds)[,sample_col] == y])))
    }
    return(tmp)
  }
  get_sample_means <- function(sample_names) {
    tmp <- as.data.frame(lapply(sample_names, function(sample) {
      tmp1 <- c(unlist(lapply(group_values, get_group_mean, 
                              cds = cds, y = sample)))
    }))
    sample_means <- as.data.frame(cbind(group_values, tmp))
    colnames(sample_means) <- c("group", sprintf("%s_means", 
                                                 samples))
    return(sample_means)
  }
  
  sample_means <- get_sample_means(samples)
  
  pData(cds)$sample_mean <- 0
  for (i in samples) {
    pData(cds)[pData(cds)[sample_col] == i, ]$sample_mean <- as.numeric(as.character(sample_means[match(pData(cds)[pData(cds)[,sample_col] == i, ][,group], sample_means$group), paste(i, "_means",  sep = "")]))
  }
  sample_means$Overall_mean <- lapply(group_values, get_group_mean, 
                                      cds = cds)
  
  pData(cds)$newCDS_mean <- as.numeric(as.character(sample_means$Overall_mean[match(pData(cds)[,group], 
                                                                                    sample_means$group)]))
  
  newCDS_exprs <- Matrix::t((Matrix::t(exprs(cds))/pData(cds)$sample_mean) * pData(cds)$newCDS_mean)
  
  newCDS <- newCellDataSet(as(as.matrix(newCDS_exprs), "sparseMatrix"), 
                          phenoData = new("AnnotatedDataFrame", data = pData(cds)), 
                         featureData = new("AnnotatedDataFrame", data = fData(cds)), 
                        lowerDetectionLimit = 1, expressionFamily = VGAM::negbinomial.size())
  return(newCDS_exprs)
}