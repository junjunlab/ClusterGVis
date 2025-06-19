#' Check and Install Required Packages
#'
#' This function checks if the required Bioconductor and CRAN packages are installed.
#' If any of the required packages are not installed, it installs them automatically.
#' This includes Bioconductor packages such as \code{clusterProfiler}, \code{SummarizedExperiment},
#' \code{ComplexHeatmap}, \code{SingleCellExperiment}, \code{TCseq}, \code{monocle}, and \code{Biobase}.
#' It also checks for CRAN packages such as \code{Seurat}, \code{WGCNA}, \code{igraph},
#' \code{pheatmap}, \code{circlize}, \code{e1071}, and \code{colorRamps}.
#'
#' @return NULL This function does not return a value. It installs the packages as needed.
#' @importFrom BiocManager install
#' @importFrom utils install.packages
#' @export
check_dependency <- function(){
  # installing needed packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager")
  }

  packages <- c("clusterProfiler",
                "SummarizedExperiment",
                "ComplexHeatmap",
                "SingleCellExperiment",
                "TCseq","monocle","Biobase")

  # check needed packages
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }

  crans <- c("Seurat","WGCNA","igraph","pheatmap","circlize","e1071","colorRamps")

  # check needed packages
  for (pkg in crans) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      utils::install.packages(pkg)
    }
  }
}
