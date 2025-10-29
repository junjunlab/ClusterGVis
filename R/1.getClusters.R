#' Determine Optimal Clusters for Gene Expression or Pseudotime Data
#'
#' @author JunZhang
#'
#' The `getClusters` function identifies the optimal number of clusters for a
#' given data object.
#' It supports multiple input types, including gene expression matrices and
#' objects such as `cell_data_set`. The function implements the Elbow method
#' to evaluate within-cluster sum of squares (WSS) across a range of cluster
#' numbers and visualizes the results.
#'
#' @param obj A data object representing the gene expression data or
#' pseudotime data:
#'   - If the input is a `cell_data_set` object (e.g., from `Monocle3`),
#'     the function preprocesses the data using `pre_pseudotime_matrix`.
#'   - If the input is a numeric matrix or a `data.frame`, it directly uses
#'   this data, or SummarizedExperiment object.
#'   Default is `NULL`.
#' @param ... Additional arguments passed to the preprocessing function
#'   `pre_pseudotime_matrix` (e.g., `assays`, `normalize`, etc.).
#'
#' @return A `ggplot` object visualizing the Elbow plot, where:
#'   - The x-axis represents the number of clusters tested.
#'   - The y-axis represents the WSS for each cluster number.
#'
#' The optimal cluster number can be visually identified at the "elbow point,
#' " where the
#' reduction in WSS diminishes sharply.
#'
#' @examples
#'
#' data("exps")
#' getClusters(obj = exps)
#'
#' @return a ggplot.
#' @export
#'
getClusters <- function(obj = NULL, ...) {
  # check datatype
  cls <- class(obj)

  if ("cell_data_set" %in% cls) {
    extra_params <- list(cds_obj = obj, assays = "counts", ...)
    exp <- do.call(pre_pseudotime_matrix, extra_params)
  } else if ("matrix" %in% cls | "data.frame" %in% cls) {
    exp <- obj
  } else if ("SummarizedExperiment" %in% cls){
    exp <- SummarizedExperiment::assay(obj)
  } else{
    stop("Please supply a data format with cell_data_set object,
         SummarizedExperiment object, matrix or data.frame!")
  }

  factoextra::fviz_nbclust(exp, stats::kmeans, method = "wss") +
    ggplot2::labs(subtitle = "Elbow method")
}
