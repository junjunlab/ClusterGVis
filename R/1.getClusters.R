#' Determine Optimal Clusters for Gene Expression or Pseudotime Data
#'
#' @author JunZhang
#'
#' The `getClusters` function identifies the optimal number of clusters for a given data object.
#' It supports multiple input types, including gene expression matrices and objects such as
#' `cell_data_set`. The function implements the Elbow method to evaluate within-cluster
#' sum of squares (WSS) across a range of cluster numbers and visualizes the results.
#'
#' @param obj A data object representing the gene expression data or pseudotime data:
#'   - If the input is a `cell_data_set` object (e.g., from `Monocle3`),
#'     the function preprocesses the data using `pre_pseudotime_matrix`.
#'   - If the input is a numeric matrix or a `data.frame`, it directly uses this data.
#'   Default is `NULL`.
#' @param ... Additional arguments passed to the preprocessing function
#'   `pre_pseudotime_matrix` (e.g., `assays`, `normalize`, etc.).
#'
#' @return A `ggplot` object visualizing the Elbow plot, where:
#'   - The x-axis represents the number of clusters tested.
#'   - The y-axis represents the WSS for each cluster number.
#'
#' The optimal cluster number can be visually identified at the "elbow point," where the
#' reduction in WSS diminishes sharply.
#'
#' @return a ggplot.
#' @export
#'
#' @examples
#' \donttest{
#' data("exps")
#' getClusters(exps)
#' }
getClusters <- function(obj = NULL,...){
  # check datatype
  cls <- class(obj)

  if(cls == "cell_data_set"){
    extra_params <- list(cds_obj = obj,assays = "counts",...)
    exp <- do.call(pre_pseudotime_matrix,extra_params)
  }else if(cls %in% c("matrix","data.frame")){
    exp <- obj
  }

  factoextra::fviz_nbclust(exp, stats::kmeans, method = "wss") +
    ggplot2::labs(subtitle = "Elbow method")
}
