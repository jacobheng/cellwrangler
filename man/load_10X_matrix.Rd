% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_10X_matrix.R
\name{load_10X_matrix}
\alias{load_10X_matrix}
\title{Load an expression matrix from a directory of cellranger-aligned 10X Genomics scRNAseq data}
\usage{
load_10X_matrix(cellranger_outs_path, which_matrix = "raw",
  cellranger_v3 = T)
}
\arguments{
\item{cellranger_outs_path}{the path to the "outs" directory in the cellranger library folder e.g.
"/mycellranger_library/outs"}

\item{which_matrix}{must be "raw" or "filtered" to specify the raw or filtered matrix to load respectively.}

\item{cellranger_v3}{logical; if cellranger version 3 and above was used to produced matrices.}

\item{cellranger_filter}{whether to use the filtered matrices from the cellranger pipeline as the 
expression matrix (i.e. as cells)}
}
\value{
a dgTMatrix
}
\description{
load_10X_matrix() loads an expression matrix as a sparseMatrix from a directory of 10X Genomics 
scRNAseq data created using the cellranger pipeline. Either the raw or the filtered matrix can be loaded.
}
\examples{
raw_mtx <- load_10X_matrix("/mycellranger_library/outs", which_matrix="raw")
}
\keyword{load_10X_matrix}
