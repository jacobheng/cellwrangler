% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dim_coord_barcode.R
\name{dim_coord_barcode}
\alias{dim_coord_barcode}
\title{Find rownames of cells from coordinates of dimensionality reduction}
\usage{
dim_coord_barcode(dim_proj, x1, x2, y1, y2)
}
\arguments{
\item{x1}{lower limit of x-coordinates}

\item{x2}{higher limit of x-coordinates}

\item{y1}{lower limit of y-coordinates}

\item{y2}{higher limit of y-coordinates}

\item{dim_reduction}{a dataframe specifying the coordinates of a dimensionality reduction output (e.g. 
t-SNE, UMAP etc.). Each column should specify one dimension. Rownames should correspond to some of identifier
for cells (e.g. barcode).}
}
\value{
A vector of rownames of dim_reduction
}
\description{
dim_coord_barcode() returns the rownames of cells based on the coordiates of a dimensionality reduction
projection (e.g. t-SNE, UMAP).
}
\examples{
dim_coord_barcode(UMAP_proj, 1, 1.5, 2, 3)
}
\keyword{dim_coord_barcode}
