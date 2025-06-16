#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL



# from mfuzz
standardise <- function(eset){
  data <- Biobase::exprs(eset)
  for (i in 1:dim(data)[[1]]){
    data[i,] <- (data[i,] - base::mean(data[i,],na.rm=TRUE))/stats::sd(data[i,],na.rm=TRUE)
  }
  Biobase::exprs(eset) <- data
  eset
}


# from mfuzz
mestimate<- function(eset){
  N <-  dim(Biobase::exprs(eset))[[1]]
  D <- dim(Biobase::exprs(eset))[[2]]
  m.sj <- 1 + (1418/N + 22.05)*D^(-2) + (12.33/N +0.243)*D^(-0.0406*log(N) - 0.1134)
  return(m.sj)
}


# Test whether a matrix is one of our supported sparse matrices
# author https://github.com/cole-trapnell-lab/monocle3
is_sparse_matrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix")
}

#' Get the size factors from a cds object.
#'
#' A wrapper around \code{colData(cds)$Size_Factor}
#'
#' @param cds A cell_data_set object.
#' @return An updated cell_data_set object
#'
#' @export
size_factors <- function( cds ) {
  stopifnot( methods::is( cds, "cell_data_set" ) )
  sf <- SingleCellExperiment::colData(cds)$Size_Factor
  names( sf ) <- colnames(SingleCellExperiment::counts(cds) )
  sf
}




#' The cell_data_set class from https://github.com/cole-trapnell-lab/monocle3/blob/master/R/cell_data_set.R
#'
#' The main class used by Monocle3 to hold single-cell expression data.
#' cell_data_set extends the Bioconductor SingleCellExperiment class.
#'
#' This class is initialized from a matrix of expression values along with cell
#' and feature metadata.
#'
#' @field reduce_dim_aux SimpleList, auxiliary information from reduced
#'   dimension.
#' @field principal_graph_aux SimpleList, auxiliary information from principal
#'   graph construction
#' @field principal_graph SimpleList of igraph objects containing principal
#'   graphs for different dimensionality reduction.
#' @field clusters SimpleList of cluster information for different
#'   dimensionality reduction.
#' @name cell_data_set
#' @rdname cell_data_set
#' @aliases cell_data_set-class
#' @exportClass cell_data_set
#' @importFrom Biobase package.version
#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData
#' @importFrom SingleCellExperiment reducedDim<- reducedDim reducedDims<-
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment Assays colData<- rowData<- assays assays<-
setClass("cell_data_set",
         contains = c("SingleCellExperiment"),
         slots = c(reduce_dim_aux = "SimpleList",
                   principal_graph_aux="SimpleList",
                   principal_graph = "SimpleList",
                   clusters = "SimpleList")
)


#' Generic to extract pseudotime from CDS object
#'
#' @author https://github.com/cole-trapnell-lab/monocle3
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract pseudotime for.
#'
#'
#' @return Pseudotime values.
#'
#' @export
setGeneric("pseudotime", function(x, reduction_method = "UMAP")
  standardGeneric("pseudotime"))


#' Method to extract pseudotime from CDS object
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract clusters for.
#'
#' @return Pseudotime values.
#'
#' @export
setMethod("pseudotime", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@principal_graph_aux[[
              reduction_method]]$pseudotime[colnames(exprs(x))]
            if (is.null(value)) {
              stop("No pseudotime calculated for reduction_method = ",
                   reduction_method, ". Please first run ",
                   "order_cells with reduction_method = ",
                   reduction_method, ".")
            }
            return(value)
          })


#' Generic to access cds count matrix
#'
#' @author https://github.com/cole-trapnell-lab/monocle3
#'
#' @param x A cell_data_set object.
#'
#'
#' @return Count matrix.
#'
#' @export
setGeneric("exprs", function(x) standardGeneric("exprs"))

#' Method to access cds count matrix
#' @param x A cell_data_set object.
#'
#' @return Count matrix.
#'
#' @export
setMethod("exprs", "cell_data_set", function(x) {
  value <- SummarizedExperiment::assays(x)$counts
  return(value)
})


#' Return a size-factor normalized and (optionally) log-transformed expression
#'
#' @author https://github.com/cole-trapnell-lab/monocle3
#'
#' matrix
#'
#' @param cds A CDS object to calculate normalized expression matrix from.
#' @param norm_method String indicating the normalization method. Options are
#'   "log" (Default), "binary" and "size_only".
#' @param pseudocount A pseudocount to add before log transformation. Ignored
#'   if norm_method is not "log". Default is 1.
#' @return Size-factor normalized, and optionally log-transformed, expression
#'   matrix.
#'
#'
#' @importFrom SingleCellExperiment counts
#'
#' @export
normalized_counts <- function(cds,
                              norm_method=c("log", "binary", "size_only"),
                              pseudocount=1){
  norm_method = match.arg(norm_method)
  norm_mat = SingleCellExperiment::counts(cds)
  if (norm_method == "binary"){
    # The '+ 0' coerces the matrix to type numeric. It's possible
    # to use 'as.numeric(norm_mat > 0)' but the matrix
    # attributes disappear...
    norm_mat = (norm_mat > 0) + 0
    if (is_sparse_matrix(norm_mat)){
      norm_mat = methods::as(norm_mat, "dgCMatrix")
    }
  }
  else {
    if (is_sparse_matrix(norm_mat)){
      norm_mat@x = norm_mat@x / rep.int(size_factors(cds), diff(norm_mat@p))
      if (norm_method == "log"){
        if (pseudocount == 1){
          norm_mat@x = log10(norm_mat@x + pseudocount)
        }else{
          stop("Pseudocount must equal 1 with sparse expression matrices")
        }
      }
    }else{
      norm_mat = Matrix::t(Matrix::t(norm_mat) / size_factors(cds))
      if (norm_method == "log"){
        norm_mat@x <- log10(norm_mat + pseudocount)
      }
    }
  }
  return(norm_mat)
}
