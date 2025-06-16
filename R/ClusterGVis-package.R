#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields = "Version")
  if (requireNamespace("cli", quietly = TRUE)) {
    start_up_mg <- cli::cat_boxx("Welcome to use ClusterGVis package for analysis.",
                                 col = "#347433")
  } else {
    warning("Package 'cli' is needed for this function to work.")
  }

  packageStartupMessage(start_up_mg)
  packageStartupMessage(paste("The version of ClusterGVis:",
                              pkgVersion,
                              "\nAny advice or suggestions please contact with me: 3219030654@stu.cpu.edu.cn.\n\n",
                              sep = " "))

  # installing needed packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    utils::install.packages("BiocManager")
  }

  packages <- c("clusterProfiler",
                "SummarizedExperiment",
                "ComplexHeatmap",
                "circlize",
                "SingleCellExperiment",
                "TCseq")

  # check needed packages
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    } else {
      packageStartupMessage("Checking packages needed...")
      packageStartupMessage(paste(pkg, "has been installed!"))
    }
  }
}
