% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_highbcv_genes.R
\name{find_highbcv_genes}
\alias{find_highbcv_genes}
\title{Find genes with high biological coefficent of variation}
\usage{
find_highbcv_genes(cds, annotated_fData = NULL, scaled_matrix = NULL,
  cell_threshold = 10, exprs_threshold = 0.001,
  resid_threshold = 0.2)
}
\arguments{
\item{cds}{a monocle CellDataSet object}

\item{annotated_fData}{annotated fData object produced by the cellwrangler annotate_bcv() function; if NULL,
parameters of bcv will be calculated from cds provided.}

\item{scaled_matrix}{a scaled expression matrix derived from cds provided. Exprs(cds) should be divided by 
the colsums and multiplied by 10,000.}

\item{cell_threshold}{threshold for expression in number of cells}

\item{exprs_threshold}{mean expression threshold for expression levels}

\item{resid_threshold}{threshold for residuals}
}
\description{
find_highbcv_genes() return genes with a high biological coefficent of variation (bcv) from a 
monocle CellDataSet object, subject to a number of thresholds including expression in number of cells, 
expression levels and residuals.
}
\examples{
find_highbcv_genes(dat)
}
\keyword{find_highbcv_genes}
