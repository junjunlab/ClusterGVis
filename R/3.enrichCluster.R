globalVariables(c("Description", "group", "pvalue", "geneID"))

#' Perform GO/KEGG Enrichment Analysis for Multiple Clusters
#' @author JunZhang
#'
#' This function performs Gene Ontology (GO) or KEGG enrichment analysis, or
#'  custom gene set enrichment,
#' on clustered genes. It supports multiple clusters, incorporating
#' cluster-specific results into its analysis.
#'
#' @param object An object containing clustering results.
#'   This is clusterData object. Alternatively, it can be a `CellDataSet`
#'   object, in which case the function can also visualize pseudotime data.
#' @param type Character. The type of enrichment analysis to perform. Options
#' include:
#'   - `"BP"`: Biological Process (GO)
#'   - `"MF"`: Molecular Function (GO)
#'   - `"CC"`: Cellular Component (GO)
#'   - `"KEGG"`: KEGG Pathway analysis
#'   - `"ownSet"`: Custom gene set enrichment, requiring `TERM2GENE` and
#'   optionally `TERM2NAME`.
#' @param TERM2GENE A data frame containing mappings of terms to genes.
#' Required when `type = "ownSet"`.
#'   This must be a two-column data frame, where the first column is the
#'   term and the second column is the gene.
#' @param TERM2NAME A data frame containing term-to-name mappings. Optional
#' when `type = "ownSet"`.
#'   This must also be a two-column data frame, where the first column is the
#'   term and the second column is the name.
#' @param OrgDb An organism database object (e.g., `org.Hs.eg.db` for human or
#' `org.Mm.eg.db` for mouse),
#'   used for GO or KEGG enrichment analysis.
#' @param id.trans Logical. Whether to perform gene ID transformation. Default
#'  is `TRUE`.
#' @param fromType Character. The type of the input gene IDs (e.g.,
#' `"SYMBOL"`, `"ENSEMBL"`). Default is `"SYMBOL"`.
#' @param toType Character. The target ID type for transformation using
#' `clusterProfiler::bitr`
#'   (e.g., `"ENTREZID"`). Default is `"ENTREZID"`.
#' @param readable Logical. Whether to convert the enrichment result IDs back
#'  to a readable format (e.g., SYMBOL).
#'   Only applicable for GO and KEGG analysis. Default is `TRUE`.
#' @param organism Character. The KEGG organism code (e.g., `"hsa"` for human,
#'  `"mmu"` for mouse). Required when
#'   performing KEGG enrichment. Default is `"hsa"`.
#' @param pvalueCutoff Numeric. The p-value cutoff for enriched terms to be
#' included in the results. Default is `0.05`.
#' @param topn Integer or vector. The number of top enrichment results to
#' extract. If a single value, it is applied
#'   to all clusters. Otherwise, it should match the number of clusters.
#'    Default is `5`.
#' @param seed Numeric. Seed for random operations to ensure reproducibility.
#'  Default is `5201314`.
#' @param add.gene Logical. Whether to include the list of genes associated
#'  with each enriched term in the results.
#'   Default is `FALSE`.
#' @param use_internal_data Logical, use KEGG.db or latest online KEGG data for
#'  enrichKEGG function.
#' Default is `FALSE`.
#' @param heatmap.type Character. The type of heatmap visualization to use when
#' input data is a `CellDataSet` object.
#'   Options include:
#'   - `"plot_pseudotime_heatmap2"`
#'   - `"plot_genes_branched_heatmap2"`
#'   - `"plot_multiple_branches_heatmap2"`
#' @param ... Additional arguments passed to plot_pseudotime_heatmap2/
#' plot_genes_branched_heatmap2/plot_multiple_branches_heatmap2 functions.
#'
#' @examples
#'
#' library(org.Mm.eg.db)
#'
#' data("exps")
#'
#' # kmeans
#' ck <- clusterData(
#'   obj = exps,
#'   cluster.method = "kmeans",
#'   cluster.num = 3
#' )
#'
#' # enrich for clusters
#' enrich <- enrichCluster(
#'   object = ck,
#'   OrgDb = org.Mm.eg.db,
#'   type = "BP",
#'   pvalueCutoff = 0.05,
#'   topn = 5
#' )
#'
#' # check
#' head(enrich, 3)
#'
#' @return a data.frame.
#' @export
#'
enrichCluster <- function(object = NULL,
                          type = c("BP", "MF", "CC", "KEGG", "ownSet"),
                          TERM2GENE = NULL,
                          TERM2NAME = NULL,
                          OrgDb = NULL,
                          id.trans = TRUE,
                          fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          readable = TRUE,
                          # for kegg mouse(mmu)
                          organism = "hsa",
                          pvalueCutoff = 0.05,
                          topn = 5,
                          seed = 5201314,
                          add.gene = FALSE,
                          use_internal_data = FALSE,
                          heatmap.type = c(
                            "plot_pseudotime_heatmap2",
                            "plot_genes_branched_heatmap2",
                            "plot_multiple_branches_heatmap2"
                          ),
                          ...) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required. Please install it.")
  }

  type <- match.arg(type, c("BP", "MF", "CC", "KEGG", "ownSet"))
  heatmap.type <- match.arg(
    heatmap.type,
    c(
      "plot_pseudotime_heatmap2",
      "plot_genes_branched_heatmap2",
      "plot_multiple_branches_heatmap2"
    )
  )

  # check datatype
  cls <- class(object)

  # check heatmap.type
  if (cls == "CellDataSet") {
    if (heatmap.type %in% c(
      "plot_pseudotime_heatmap2",
      "plot_genes_branched_heatmap2"
    )) {
      extra_params <- list(cds_subset = object, ...)
      object <- do.call(plot_pseudotime_heatmap2, extra_params)
    } else {
      extra_params <- list(cds = object, ...)
      object <- do.call(plot_pseudotime_heatmap2, extra_params)
    }
  }

  # get data
  enrich.data <- object$wide.res

  # check geneType
  if (startsWith(object$geneType, "nounique")) {
    split <- unlist(strsplit(object$geneType, split = "\\|"))[2]
    # enrich.data$gene <- sapply(strsplit(as.character(enrich.data$gene),
    # split = split),"[",1)
    enrich.data$gene <- vapply(strsplit(as.character(enrich.data$gene),
                                        split = split), function(x) {
      x[1]
    }, character(1))
  }

  # loop for enrich
  purrr::map_df(seq_len(length(unique(
    enrich.data$cluster
  ))), function(x) {
    # filter
    tmp <- enrich.data |> 
      dplyr::filter(cluster == unique(enrich.data$cluster)[x])

    # =============================================
    # enrich

    # entriz id transformation
    # id.trans <- ifelse(type == "ownSet",FALSE,TRUE)

    if (id.trans == TRUE) {
      gene.ent <- clusterProfiler::bitr(
        tmp$gene,
        fromType = fromType,
        toType = toType,
        OrgDb = OrgDb
      )

      tartget.gene <- unlist(gene.ent[, c(toType)])
    } else {
      tartget.gene <- tmp$gene
    }

    # GO enrich
    if (type %in% c("BP", "MF", "CC")) {
      # set.seed(seed)
      ego <- clusterProfiler::enrichGO(
        gene          = tartget.gene,
        keyType       = toType,
        OrgDb         = OrgDb,
        ont           = type,
        pAdjustMethod = "BH",
        pvalueCutoff  = 1,
        qvalueCutoff  = 1,
        readable      = readable
      )
    } else if (type == "ownSet") {
      # set.seed(seed)
      ego <- clusterProfiler::enricher(
        gene = tartget.gene,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1
      )
    } else {
      # set.seed(seed)
      ego <- clusterProfiler::enrichKEGG(
        gene = tartget.gene,
        keyType = "kegg",
        organism = organism,
        universe = NULL,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        use_internal_data = use_internal_data
      )

      # transform gene id
      if (readable == TRUE) {
        ego <- clusterProfiler::setReadable(ego, OrgDb = OrgDb,
                                            keyType = toType)
      } else {
        ego <- ego
      }
    }

    # to data.frame
    enrich_res <- data.frame(ego)
    if (nrow(enrich_res) > 0) {
      df <- data.frame(ego) |> 
        dplyr::filter(pvalue < pvalueCutoff) |> 
        dplyr::mutate(group = paste("C", unique(enrich.data$cluster)[x],
                                    sep = "")) |> 
        dplyr::arrange(pvalue)

      # whether save all res
      if (!is.null(topn)) {
        # top n selection
        if (length(topn) == 1) {
          top <- topn
        } else {
          top <- topn[x]
        }

        df <- df |> 
          # dplyr::select(group,Description,pvalue) |> 
          dplyr::slice_head(n = top)

        # add gene enrich ratio
        df1 <- purrr::map_df(seq_len(nrow(df)), function(x) {
          tmp1 <- df[x, ]
          size <- unlist(strsplit(as.character(tmp1$GeneRatio), split = "/"))
          tmp1$ratio <- (as.numeric(size[1]) / as.numeric(size[2])) * 100
          return(tmp1)
        })

        # select columns
        if (add.gene == TRUE) {
          df2 <- df1 |>dplyr::select(group, Description, pvalue,
                                       ratio, geneID)
        } else {
          df2 <- df1 |>dplyr::select(group, Description, pvalue, ratio)
        }
      } else {
        df2 <- df
      }

      # remove no pass threshold terms
      df2 <- df2 |>stats::na.omit()

      # results
      return(df2)
    } else {
      return(NULL)
    }
  })
}
