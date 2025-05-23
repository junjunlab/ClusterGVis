% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pseudotime}
\alias{pseudotime}
\title{Generic to extract pseudotime from CDS object}
\usage{
pseudotime(x, reduction_method = "UMAP")
}
\arguments{
\item{x}{A cell_data_set object.}

\item{reduction_method}{Reduced dimension to extract pseudotime for.}
}
\value{
Pseudotime values.
}
\description{
Generic to extract pseudotime from CDS object
}
\author{
https://github.com/cole-trapnell-lab/monocle3
}
