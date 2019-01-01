#'Find rownames of cells from coordinates of dimensionality reduction
#'
#' dim_coord_barcode() returns the rownames of cells based on the coordiates of a dimensionality reduction
#' projection (e.g. t-SNE, UMAP).
#' 
#' @param dim_reduction a dataframe specifying the coordinates of a dimensionality reduction method (e.g. 
#' t-SNE, UMAP etc.). Each column should specify one dimension. Rownames should correspond to some of identifier
#' for cells (e.g. barcode).
#' @param x1 lower limit of x-coordinates
#' @param x2 higher limit of x-coordinates
#' @param y1 lower limit of y-coordinates
#' @param y2 higher limit of y-coordinates
#' @keywords dim_coord_barcode
#' @export
#' @return A vector of rownames of dim_reduction
#' @examples
#' dim_coord_barcode(UMAP_proj, 1, 1.5, 2, 3)

#dim_coord_barcode function
dim_coord_barcode <- function(dim_proj,x1,x2,y1,y2) {
  rownames(dim_proj[dim_proj[, 1] >= x1 & dim_proj[, 1] <= x2 & dim_proj[, 2] >= y1 & dim_proj[, 2] <= y2, ])
}