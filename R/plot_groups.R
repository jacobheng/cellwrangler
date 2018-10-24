#'Plot cell groups
#'
#' plot_groups() colors cells differently according to a grouping variable (e.g. cell clusters) on 
#' a two-dimensional plot. A two-dimensional projection object with coordinates for each cell (e.g. by t-SNE or UMAP dimensionality 
#' reduction) is required. 
#' @param group_vector a vector assigning each cell to a group
#' @param projection a dataframe or matrix containing x and y coordinates for each cell in two columns
#' (e.g. derived from t-SNE or UMAP)
#' @param colors a vector of colors to use for coloring cells. Number of colors must equal to the number 
#' of clusters specified in cluster_vector
#' @param rescale logical; whether to rescale color scale as a continuous gradient from the minimum to 
#' the maximum value of expression for specified gene; defaults to FALSE
#' @param cell_size, the size of the point representing each cell
#' @keywords plot_groups
#' @export
#' @return A ggplot2 object
#' @examples
#' plot_groups(group_vector=kmeans_clusters_res, projection=tsne_proj[,c("TSNE.1", "TSNE.2")])


plot_groups <- function (group_vector, projection, colors = NULL, alpha = 1,
                               cell_size = 0.1, title = NULL)
{
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_group <- data.frame(cbind(projection, group_vector))
  proj_cols <- colnames(projection)
  proj_group[proj_cols] <- lapply(proj_group[proj_cols], function(x){
    num.x <- as.numeric(as.character(x))
    return(num.x)
  } )
  proj_group_melt <- melt(proj_group, id.vars = c("Component.1",
                                              "Component.2"))
  print(head(proj_group_melt))
  colnames(proj_group_melt) <-  c("Component.1","Component.2", "group_vector", "group")
  p <- ggplot(proj_group_melt, aes(Component.1, Component.2)) +
    geom_point(aes(colour = value), size = cell_size,
               alpha = alpha) + guides(col = guide_legend(title = "Group",
                                                          override.aes = list(size = 3))) + facet_wrap(~variable) +
    labs(x = projection_names[1], y = projection_names[2]) +
    theme_bw()
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if ((!is.null(colors))) {
    names(colors) <- 1:length(colours)
    p <- p + scale_color_manual(values = colors)
  }
  p <- p + theme(plot.title = element_text(hjust = 0.5),
                 strip.text.x = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
  return(p)
}

