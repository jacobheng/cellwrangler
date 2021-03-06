% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_pca_biplot.R
\name{sparse_pca_biplot}
\alias{sparse_pca_biplot}
\title{Plot biplot for sparse pca}
\usage{
sparse_pca_biplot(sparse_pca_obj, pcs, group = NULL)
}
\arguments{
\item{sparse_pca_obj}{a sparse pca object created by the cellRangerRkit sparse_pca() function}

\item{pcs}{a two-component vector of PCs to plot e.g. c(1,2)}

\item{group}{a vector specifying groups for the observations. The points on the biplot will be colored
according to group}
}
\value{
A ggplot object
}
\description{
sparse_pca_biplot() creates a biplot for two principal components computed by sparse PCA 
(implemented by cellrangerRkit).
}
\examples{
sparse_pca_biplot(sparse_pca_obj, pcs=c(1,2), group=pData(dat)$CellType)
}
\keyword{sparse_pca_biplot()}
