#' @name getClusters
#' @author JunZhang
#' @title using getClusters to define a optimal cluster numbers
#'
#' @param exp the gene normalized expression data(fpkm/rpkm/tpm), default NULL.
#'
#' @return a ggplot.
#' @export
#'
#' @examples
#' \dontrun{
#' data("exps")
#' getClusters(exps)
#' }
getClusters <- function(exp = NULL){
  factoextra::fviz_nbclust(exp, stats::kmeans, method = "wss") +
    ggplot2::labs(subtitle = "Elbow method")
}
