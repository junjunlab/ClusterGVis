


# from mfuzz
standardise <- function(eset) {
  data <- Biobase::exprs(eset)
  for (i in seq_len(dim(data)[[1]])) {
    data[i, ] <-
      (data[i, ] - base::mean(data[i, ], na.rm = TRUE)) / stats::sd(data[i, ],
        na.rm =
          TRUE
      )
  }
  Biobase::exprs(eset) <- data
  eset
}


# from mfuzz
mestimate <- function(eset) {
  N <- dim(Biobase::exprs(eset))[[1]]
  D <- dim(Biobase::exprs(eset))[[2]]
  m.sj <-
    1 + (1418 / N + 22.05) * D^(-2) + (12.33 / N + 0.243) * D^(-0.0406 *
      log(N) - 0.1134)
  return(m.sj)
}


# Test whether a matrix is one of our supported sparse matrices
# author https://github.com/cole-trapnell-lab/monocle3
is_sparse_matrix <- function(x) {
  class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix")
}



size_factors <- function(cds) {
  stopifnot(methods::is(cds, "cell_data_set"))

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required. Please install it.")
  }

  sf <- SingleCellExperiment::colData(cds)$Size_Factor
  names(sf) <- colnames(SingleCellExperiment::counts(cds))
  sf
}




# ==============================================================================

#' Cell Data Set Class
#' @name cell_data_set-class
#' @rdname cell_data_set-class
#' @exportClass cell_data_set
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setClass(
  "cell_data_set",
  contains = "SingleCellExperiment",
  slots = c(
    reduce_dim_aux = "SimpleList",
    principal_graph_aux = "SimpleList",
    principal_graph = "SimpleList",
    clusters = "SimpleList"
  )
)



# ==============================================================================

#' Generic to extract pseudotime from CDS object
#'
#' @author https://github.com/cole-trapnell-lab/monocle3
#'
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract pseudotime for.
#'
#'
#' @return A named numeric vector containing pseudotime values for each cell.
#'   The vector has the following properties:
#'   \itemize{
#'     \item Names correspond to cell identifiers (column names from the
#'     expression matrix)
#'     \item Values are numeric pseudotime values representing trajectory
#'     progress
#'     \item Values typically range from 0 (trajectory start) to maximum
#'     trajectory length
#'     \item Higher values indicate cells further along the developmental
#'      trajectory
#'   }
#'
#'
#'
setGeneric("pseudotime", function(x, reduction_method = "UMAP") {
  standardGeneric("pseudotime")
})



#' @rdname pseudotime
#'
setMethod("pseudotime", "cell_data_set",
          function(x, reduction_method = "UMAP") {
  value <- x@principal_graph_aux[[reduction_method]]$pseudotime[
    colnames(exprs(x))]
  if (is.null(value)) {
    stop(
      "No pseudotime calculated for reduction_method = ",
      reduction_method,
      ". Please first run ",
      "order_cells with reduction_method = ",
      reduction_method,
      "."
    )
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
#'
setGeneric("exprs", function(x) {
  standardGeneric("exprs")
})

#' Method to access cds count matrix
#' @param x A cell_data_set object.
#'
#' @return Count matrix.
#'
#'
setMethod("exprs", "cell_data_set", function(x) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required. Please install it.")
  }

  value <- SummarizedExperiment::assays(x)$counts
  return(value)
})


# ==============================================================================

# Return a size-factor normalized and (optionally) log-transformed expression
#
# @author https://github.com/cole-trapnell-lab/monocle3
#
normalized_counts <- function(cds,
                              norm_method = c("log", "binary", "size_only"),
                              pseudocount = 1) {
  norm_method <- match.arg(norm_method,
                           choices = c("log", "binary", "size_only"))

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required. Please install it.")
  }

  norm_mat <- SingleCellExperiment::counts(cds)
  if (norm_method == "binary") {
    # The '+ 0' coerces the matrix to type numeric. It's possible
    # to use 'as.numeric(norm_mat > 0)' but the matrix
    # attributes disappear...
    norm_mat <- (norm_mat > 0) + 0
    if (is_sparse_matrix(norm_mat)) {
      norm_mat <- methods::as(norm_mat, "dgCMatrix")
    }
  } else {
    if (is_sparse_matrix(norm_mat)) {
      norm_mat@x <- norm_mat@x / rep.int(size_factors(cds), diff(norm_mat@p))
      if (norm_method == "log") {
        if (pseudocount == 1) {
          norm_mat@x <- log10(norm_mat@x + pseudocount)
        } else {
          stop("Pseudocount must equal 1 with sparse expression matrices")
        }
      }
    } else {
      norm_mat <- Matrix::t(Matrix::t(norm_mat) / size_factors(cds))
      if (norm_method == "log") {
        norm_mat@x <- log10(norm_mat + pseudocount)
      }
    }
  }
  return(norm_mat)
}
