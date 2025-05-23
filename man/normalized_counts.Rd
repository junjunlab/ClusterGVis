% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{normalized_counts}
\alias{normalized_counts}
\title{Return a size-factor normalized and (optionally) log-transformed expression}
\usage{
normalized_counts(
  cds,
  norm_method = c("log", "binary", "size_only"),
  pseudocount = 1
)
}
\arguments{
\item{cds}{A CDS object to calculate normalized expression matrix from.}

\item{norm_method}{String indicating the normalization method. Options are
"log" (Default), "binary" and "size_only".}

\item{pseudocount}{A pseudocount to add before log transformation. Ignored
if norm_method is not "log". Default is 1.}
}
\value{
Size-factor normalized, and optionally log-transformed, expression
matrix.
}
\description{
Return a size-factor normalized and (optionally) log-transformed expression
}
\author{
https://github.com/cole-trapnell-lab/monocle3

matrix
}
