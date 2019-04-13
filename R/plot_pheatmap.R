#' Plot pheatmap object
#'
#' @description Function for plotting a heatmap from an object created by assignment of output from the pheatmap
#' function in the popular R package pheatmap
#' @param pheatmap_object object containing output from the pheatmap function. Created by assigning output of 
#' pheatmap function to an object (i.e. pheatmap_object <- pheatmap(matrix))
#' @keywords plot_pheatmap
#' @export
#' @return A heatmap created by pheatmap
#' @examples
#' pheatmap_object <- pheatmap(matrix)
#' plot_pheatmap(pheatmap_object)

plot_pheatmap<-function(pheatmap_object){
  grid::grid.newpage()
  grid::grid.draw(pheatmap_object$gtable)
}