#' Prepare scRNA Data for clusterGvis Analysis
#'
#' This function prepares single-cell RNA sequencing (scRNA-seq) data for
#' differential
#' gene expression analysis. It extracts the expression data for the specified
#' cells
#' and genes, and organizes them into a dataframe format suitable for downstream
#'  analysis.
#'
#' @param object an object of class Seurat containing the scRNA-seq data.
#' @param diffData a dataframe containing information about the differential
#' expression analysis which can
#' be output from function FindAllMarkers.
#' @param showAverage a logical indicating whether to show the average gene
#' expression across all cells.
#' @param cells a vector of cell names to extract from the Seurat object. If
#' NULL, all cells will be used.
#' @param group.by a string specifying the grouping variable for differential
#' expression analysis. Default is 'ident', which groups cells by their
#' assigned clusters.
#' @param assays a string or vector of strings specifying the assay(s) to
#' extract from the Seurat object. Default is 'RNA'.
#' @param slot a string specifying the slot name where the assay data is stored
#' in the Seurat object. Default is 'data'.
#' @param scale.data whether do Z-score for expression data, default TRUE.
#' @param cluster.order the celltype orders.
#' @param keep.uniqGene a logical indicating whether to keep only unique gene
#' names. Default is TRUE.
#' @param sep a character string to separate gene and cell names in the output
#' dataframe. Default is "_".
#'
#' @return a dataframe containing the expression data for the specified genes
#' and cells,
#' organized in a format suitable for differential gene expression analysis.
#'
#' @examples
#'
#' # Load required libraries
#' library(Seurat)
#' library(dplyr)
#'
#' data("pbmc_subset")
#'
#' # report only the positive ones
#' pbmc.markers.all <- Seurat::FindAllMarkers(pbmc_subset,
#'   only.pos = TRUE,
#'   min.pct = 0.25,
#'   logfc.threshold = 0.25
#' )
#'
#' # get top 10 genes
#' pbmc.markers <- pbmc.markers.all |>
#'   dplyr::group_by(cluster) |>
#'   dplyr::top_n(n = 20, wt = avg_log2FC)
#'
#'
#' # prepare data from seurat object
#' st.data <- prepareDataFromscRNA(
#'   object = pbmc_subset,
#'   diffData = pbmc.markers,
#'   showAverage = TRUE
#' )
#'
#' @export
prepareDataFromscRNA <- function(object = NULL,
                                 diffData = NULL,
                                 showAverage = TRUE,
                                 cells = NULL,
                                 group.by = "ident",
                                 assays = "RNA",
                                 slot = "data",
                                 scale.data = TRUE,
                                 cluster.order = NULL,
                                 keep.uniqGene = TRUE,
                                 sep = "_") {
  # ============================================================================
  # get data form object
  # ============================================================================
  markerGene <- unique(diffData$gene)

  # choose mode
  if (showAverage == TRUE) {
    # get cells mean gene expression
    vr <- utils::compareVersion(
      as.character(utils::packageVersion("Seurat")), "5")
    if (vr == 1) {
      mean_gene_exp <- Seurat::AverageExpression(
        object,
        features = markerGene,
        group.by = group.by,
        assays = assays,
        layer = slot
      ) |>
        data.frame() |>
        as.matrix()
    } else {
      mean_gene_exp <- Seurat::AverageExpression(
        object,
        features = markerGene,
        group.by = group.by,
        assays = assays,
        slot = slot
      ) |>
        data.frame() |>
        as.matrix()
    }


    # add colnames
    name1 <- gsub(
      pattern = paste0(assays, ".", sep = ""),
      replacement = "",
      colnames(mean_gene_exp)
    )
    colnames(mean_gene_exp) <- gsub(pattern = "\\.", replacement = " ", name1)

    # assign colnames
    colnames(mean_gene_exp) <- levels(Seurat::Idents(object))

    # whether do zscore
    if (scale.data == TRUE) {
      mean_gene_exp <- t(scale(t(mean_gene_exp)))
    }

    # cell type orders
    if (!is.null(cluster.order)) {
      mean_gene_exp <- mean_gene_exp[, cluster.order]
    }

    geneMode <- "average"
  } else {
    # cell inro
    cell.order <- data.frame(
      cell.id = names(Seurat::Idents(object)),
      cell.ident = Seurat::Idents(object)
    )

    # order cell type
    if (is.null(cluster.order)) {
      cell.order$cell.ident <- factor(cell.order$cell.ident,
                                      levels = levels(Seurat::Idents(object)))
    } else {
      cell.order$cell.ident <- factor(cell.order$cell.ident,
                                      levels = cluster.order)
    }
    cell.order <- cell.order[order(cell.order$cell.ident), ]

    # get all cells data
    getassy <- Seurat::GetAssayData(object = object, slot = slot)[
      features = markerGene, cells = NULL, drop = FALSE] |>
      as.matrix()

    # reorder cells
    id.order <- match(cell.order$cell.id, colnames(getassy))
    getassy <- getassy[, id.order]

    # re-assign colnames
    colnames(getassy) <- paste(colnames(getassy),
                               cell.order$cell.ident, sep = "|")

    mean_gene_exp <- getassy

    # whether do zscore
    if (scale.data == TRUE) {
      mean_gene_exp <- t(scale(t(mean_gene_exp)))
    }

    geneMode <- "all"
  }

  # ============================================================================
  # prepare data
  # ============================================================================
  # add gene column
  merMat <- data.frame(mean_gene_exp, check.names = FALSE)
  merMat$gene <- rownames(merMat)
   

  # count marker gene numbers for each cluster
  cinfo.gene <- diffData[, c("cluster", "gene")]

  # loop
  cn <- unique(cinfo.gene$cluster)
  purrr::map_df(seq_along(cn), function(x) {
    tmp <- cinfo.gene[which(cinfo.gene$cluster == cn[x]), ]

    # filter data
    tmp2 <- merMat[which(merMat$gene %in% tmp$gene), ] |>
      dplyr::mutate(cluster = as.character(x))

    return(tmp2)
  }) -> wide.res

  # whether retain unique gene name
  if (keep.uniqGene == TRUE) {
    # wide.res <- wide.res |> dplyr::distinct(., gene, .keep_all = TRUE)
    # geneType <- paste("unique", sep, sep = "|")
    
    duplicated_genes <- duplicated(wide.res$gene)
    wide.res <- wide.res[!duplicated_genes, ]
    geneType <- paste0("unique", "|", sep)
  } else {
    # wide.res <- wide.res |> dplyr::mutate(.,
    #     gene = make.unique(gene, sep = sep))
    # geneType <- paste("nounique", sep, sep = "|")
    
    wide.res$gene <- make.unique(wide.res$gene, sep = sep)
    geneType <- paste0("nounique", "|", sep)
  }

  # wide to long
  df <- reshape2::melt(
    wide.res,
    id.vars = c("cluster", "gene"),
    variable.name = "cell_type",
    value.name = "norm_value"
  )

  # add cluster name
  df$cluster_name <- paste("cluster ", df$cluster, sep = "")

  if (showAverage == FALSE) {
    # df$cell_type <- sapply(strsplit(as.character(df$cell_type),
    # split = "\\|"),"[",2)
    df$cell_type <- vapply(strsplit(as.character(df$cell_type), split = "\\|"),
                           function(x) {
      x[2]
    }, character(1))
  }

  # add gene number
  # cltn <- table(wide.res$cluster)
  cl.info <- data.frame(table(wide.res$cluster)) |>
    dplyr::mutate(Var1 = as.numeric(as.character(Var1))) |>
    dplyr::arrange(Var1)

  id <- unique(df$cluster_name)
  purrr::map_df(seq_along(id), function(x) {
    tmp <- df |>
      dplyr::filter(cluster_name == id[x])

    tmp |>
      dplyr::mutate(cluster_name = paste(cluster_name,
                                         " (", cl.info$Freq[x], ")", sep = ""))
  }) -> df

  # cluster order
  df$cluster_name <- factor(df$cluster_name,
    levels = paste("cluster ", cl.info$Var1, " (", cl.info$Freq, ")", sep = "")
  )

  # return
  return(
    list(
      wide.res = wide.res,
      long.res = df,
      type = "scRNAdata",
      geneMode = geneMode,
      geneType = geneType
    )
  )
}
