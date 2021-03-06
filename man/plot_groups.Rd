% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_groups.R
\name{plot_groups}
\alias{plot_groups}
\title{Plot cell groups}
\usage{
plot_groups(group_vector, dim_reduction, colors = NULL, alpha = 1,
  cell_size = 0.1, title = NULL)
}
\arguments{
\item{group_vector}{a vector assigning each cell to a group}

\item{dim_reduction}{a dataframe specifying the coordinates of a dimensionality reduction output (e.g. 
t-SNE, UMAP etc.). Each column should specify one dimension.}

\item{colors}{a vector of colors to use for coloring cells. Number of colors must equal to the number 
of clusters specified in cluster_vector}

\item{cell_size, }{the size of the point representing each cell}

\item{rescale}{logical; whether to rescale color scale as a continuous gradient from the minimum to 
the maximum value of expression for specified gene; defaults to FALSE}
}
\value{
A ggplot2 object
}
\description{
plot_groups() colors cells differently according to a grouping variable (e.g. cell clusters) on 
a two-dimensional plot. A two-dimensional dimensionality reduction object with coordinates for each cell 
(e.g. by t-SNE or UMAP dimensionality reduction) is required.
}
\examples{
plot_groups(group_vector=kmeans_clusters_res, dim_reduction=tsne_proj[,c("TSNE.1", "TSNE.2")])
}
\keyword{plot_groups}
