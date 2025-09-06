globalVariables(
  c(
    "cell_type",
    "cluster.num",
    "gene",
    "ratio",
    "bary",
    "membership",
    "norm_value",
    "id",
    "log10P",
    "pval",
    "Var1",
    "seurat_clusters",
    "cell.ident",
    "getassy",
    "geneType",
    ".data"
  )
)

#' Visualize Clustered Gene Data Using Line Plots and Heatmaps
#'
#' This function visualizes clustered gene expression data as line plots,
#' heatmaps, or
#' a combination of both, using the `ComplexHeatmap` and `ggplot2` frameworks.
#'  Gene
#' annotations, sample annotations, and additional features like custom color
#'  schemes
#' and annotations for GO/KEGG terms are supported for visualization.
#'
#' @author JunZhang
#'
#' @title using visCluster to visualize cluster results from clusterData and
#' enrichCluster output
#'
#' @param object clusterData object, default NULL.
#' @param ht.col.list list of heatmap col_range and col_color, default
#' list(col_range = c(-2, 0, 2),col_color = c("#08519C", "white", "#A50F15")).
#' @param border whether add border for heatmap, default TRUE.
#' @param plot.type the plot type to choose which incuding "line","heatmap"
#' and "both".
#' @param ms.col membership line color form Mfuzz cluster method results,
#' default c('#0099CC','grey90','#CC3333').
#' @param line.size line size for line plot, default 0.1.
#' @param line.col line color for line plot, default "grey90".
#' @param add.mline whether add median line on plot, default TRUE.
#' @param mline.size median line size, default 2.
#' @param mline.col median line color, default "#CC3333".
#' @param ncol the columns for facet plot with line plot, default 4.
#' @param ctAnno.col the heatmap cluster annotation bar colors, default NULL.
#' @param set.md the represent line method on heatmap-line plot(mean/median),
#'  default "median".
#' @param textbox.pos the relative position of text in left-line plot,
#' default c(0.5,0.8).
#' @param textbox.size the text size of the text in left-line plot, default 8.
#' @param panel.arg the settings for the left-line panel which are
#' panel size,gap,width,fill and col, default c(2,0.25,4,"grey90",NA).
#' @param ggplot.panel.arg the settings for the ggplot2 object plot panel
#' which are
#' panel size,gap,width,fill and col, default c(2,0.25,4,"grey90",NA).
#' @param annoTerm.data the GO term annotation for the clusters, default NULL.
#' @param annoTerm.mside the wider GO term annotation box side, default "right".
#' @param termAnno.arg the settings for GO term panel annotations which
#' are fill and col,
#' default c("grey95","grey50").
#'
#' @param add.box whether add boxplot, default FALSE.
#' @param boxcol the box fill colors, default NULL.
#' @param box.arg this is related to boxplot width and border color, default
#' c(0.1,"grey50").
#' @param add.point whether add point, default FALSE.
#' @param point.arg this is related to point shape,fill,color and size,
#' default c(19,"orange","orange",1).
#' @param add.line whether add line, default TRUE.
#' @param line.side the line annotation side, default "right".
#'
#' @param markGenes the gene names to be added on plot, default NULL.
#' @param markGenes.side the gene label side, default "right".
#' @param genes.gp gene labels graphics settings, default c('italic',10,NA).
#' @param go.col the GO term text colors, default NULL.
#' @param go.size the GO term text size(numeric or "pval"), default NULL.
#' @param mulGroup to draw multiple lines annotation, supply the groups numbers
#'  with vector, default NULL.
#' @param lgd.label the lines annotation legend labels, default NULL.
#' @param show_row_names whether to show row names, default FALSE.
#' @param term.text.limit the GO term text size limit, default c(10,18).
#' @param subgroup.anno the sub-cluster for annotation, supply sub-cluster id,
#'  default NULL.
#' @param add.bar whether add bar plot for GO enrichment, default FALSE.
#' @param bar.width the GO enrichment bar width, default 8.
#' @param textbar.pos the barplot text relative position, default c(0.8,0.8).
#'
#' @param annnoblock.text whether add cluster numbers on right block
#' annotation, default TRUE.
#' @param annnoblock.gp right block annotation text color and size,
#' default c("white",8).
#' @param add.sampleanno whether add column annotation, default TRUE.
#' @param sample.group the column sample groups, default NULL.
#' @param sample.col column annotation colors, default NULL.
#' @param sample.order the orders for column samples, default NULL.
#' @param HeatmapAnnotation the 'HeatmapAnnotation' object from 'ComplexHeatmap'
#' when you have multiple annotations, default NULL.
#' @param column.split how to split the columns when supply multiple column
#' annotations, default NULL.
#' @param cluster.order the row cluster orders for user's own defination,
#' default NULL.
#' @param sample.cell.order the celltype order when input is scRNA data and
#' "showAverage = FALSE"
#' for prepareDataFromscRNA.
#' @param annoKegg.data the KEGG term annotation for the clusters, default NULL.
#' @param annoKegg.mside the wider KEGG term annotation box side, default
#'  "right".
#' @param keggAnno.arg the settings for KEGG term panel annotations which are
#'  fill and col,
#' default c("grey95","grey50").
#' @param add.kegg.bar whether add bar plot for KEGG enrichment, default FALSE.
#' @param kegg.col the KEGG term text colors, default NULL.
#' @param kegg.size the KEGG term text size(numeric or "pval"), default NULL.
#' @param by.go the GO term text box style("anno_link" or "anno_block"),
#'  default "anno_link".
#' @param by.kegg the KEGG term text box style("anno_link" or "anno_block"),
#'  default "anno_link".
#' @param word_wrap whether wrap the text, default TRUE.
#' @param add_new_line whether add new line when text is long, default TRUE.
#' @param cluster_columns whether cluster the columns, default FALSE.
#' @param pseudotime_col the branch color control for monocle input data.
#' @param gglist a list of ggplot object to annotate each cluster, default NULL.
#' @param row_annotation_obj Row annotation for heatmap, it is a
#'  `ComplexHeatmap::rowAnnotation()` object
#' when "markGenes.side" or ”line.side“ is "right". Otherwise is a list of
#' named vectors.
#'
#' @param ... othe aruguments passed by Heatmap fuction.
#'
#' @importFrom purrr map_dfr
#'
#' @return a ggplot2 or Heatmap object.
#' @export
#'
#' @examples
#'
#' data("exps")
#'
#' # mfuzz
#' cm <- clusterData(
#'   obj = exps,
#'   cluster.method = "kmeans",
#'   cluster.num = 8
#' )
#'
#' # plot
#' visCluster(
#'   object = cm,
#'   plot.type = "line"
#' )
#'
visCluster <- function(object = NULL,
                       # plot.data = NULL,
                       ht.col.list = list(
                         col_range = c(-2, 0, 2),
                         col_color = c("#08519C", "white", "#A50F15")
                       ),
                       # ht.col = c("#08519C", "white", "#A50F15"),
                       border = TRUE,
                       plot.type = c("line", "heatmap", "both"),
                       ms.col = c("#0099CC", "grey90", "#CC3333"),
                       line.size = 0.1,
                       line.col = "grey90",
                       add.mline = TRUE,
                       mline.size = 2,
                       mline.col = "#CC3333",
                       ncol = 4,
                       ctAnno.col = NULL,
                       set.md = "median",
                       textbox.pos = c(0.5, 0.8),
                       textbox.size = 8,
                       # panel size,gap,width,fill,col
                       panel.arg = c(2, 0.25, 4, "grey90", NA),
                       ggplot.panel.arg = c(2, 0.25, 4, "grey90", NA),
                       annoTerm.data = NULL,
                       annoTerm.mside = "right",
                       # textbox fill and col
                       termAnno.arg = c("grey95", "grey50"),
                       add.bar = FALSE,
                       bar.width = 8,
                       textbar.pos = c(0.8, 0.8),
                       go.col = NULL,
                       go.size = NULL,
                       by.go = "anno_link",
                       # KEGG term annotation
                       annoKegg.data = NULL,
                       annoKegg.mside = "right",
                       # textbox fill and col
                       keggAnno.arg = c("grey95", "grey50"),
                       add.kegg.bar = FALSE,
                       kegg.col = NULL,
                       kegg.size = NULL,
                       by.kegg = "anno_link",
                       word_wrap = TRUE,
                       add_new_line = TRUE,
                       # boxplot,line.point annotation
                       add.box = FALSE,
                       boxcol = NULL,
                       # box with and border color
                       box.arg = c(0.1, "grey50"),
                       add.point = FALSE,
                       # shape,fill,color,size
                       point.arg = c(19, "orange", "orange", 1),
                       add.line = TRUE,
                       line.side = "right",
                       markGenes = NULL,
                       markGenes.side = "right",
                       genes.gp = c("italic", 10, NA),
                       term.text.limit = c(10, 18),
                       mulGroup = NULL,
                       lgd.label = NULL,
                       show_row_names = FALSE,
                       subgroup.anno = NULL,
                       annnoblock.text = TRUE,
                       annnoblock.gp = c("white", 8),
                       add.sampleanno = TRUE,
                       sample.group = NULL,
                       sample.col = NULL,
                       sample.order = NULL,
                       cluster.order = NULL,
                       sample.cell.order = NULL,
                       HeatmapAnnotation = NULL,
                       column.split = NULL,
                       cluster_columns = FALSE,
                       pseudotime_col = NULL,
                       gglist = NULL,
                       row_annotation_obj = NULL,
                       ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Please install it.")
  }

  ComplexHeatmap::ht_opt(message = FALSE)

  if (is.null(ht.col.list[["col_range"]])) {
    col_range <- c(-2, 0, 2)
  } else {
    col_range <- ht.col.list[["col_range"]]
  }

  if (is.null(ht.col.list[["col_color"]])) {
    col_color <- c("#08519C", "white", "#A50F15")
  } else {
    col_color <- ht.col.list[["col_color"]]
  }

  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Please install it.")
  }

  col_fun <- circlize::colorRamp2(col_range, col_color)

  plot.type <- match.arg(plot.type)

  # choose plot type
  if (plot.type == "line") {
    # process data
    # if(is.null(plot.data)){
    #   data <- data.frame(object$long.res)
    # }else{
    #   data <- plot.data
    # }

    if (object$type %in% c("scRNAdata", "monocle", "wgcna")) {
      data <- data.frame(object$long.res)|>
        dplyr::arrange(.data[["cluster"]])
    } else {
      if (object$type %in% c("mfuzz", "TCseq")) {
        data <- data.frame(object$long.res)|>
          dplyr::arrange(.data[["cluster"]], .data[["membership"]])
      } else {
        data <- data.frame(object$long.res)|>
          dplyr::arrange(.data[["cluster"]])
      }
    }



    data$gene <- factor(data$gene, levels = unique(data$gene))

    # sample orders
    if (!is.null(sample.order)) {
      data$cell_type <- factor(data$cell_type, levels = sample.order)
    }

    # basic plot
    line <-
      ggplot2::ggplot(data, ggplot2::aes(x = cell_type, y = norm_value))

    # type
    if (object$type %in% c("mfuzz", "TCseq")) {
      line <- line +
        ggplot2::geom_line(ggplot2::aes(color = membership, group = gene),
                           size = line.size) +
        ggplot2::scale_color_gradient2(
          low = ms.col[1],
          mid = ms.col[2],
          high = ms.col[3],
          midpoint = 0.5
        )
    } else {
      line <- line +
        ggplot2::geom_line(ggplot2::aes(group = gene),
                           color = line.col,
                           size = line.size
        )
    }

    if (add.mline == TRUE) {
      if (object$type == "wgcna") {
        # line colors
        linec <- unique(data$modulecol)
        names(linec) <- linec

        line <- line +
          # median line
          ggplot2::geom_line(
            stat = "summary",
            fun = "median",
            # colour = "brown",
            size = mline.size,
            ggplot2::aes(group = 1, color = modulecol)
          ) +
          ggplot2::scale_color_manual(values = linec)
      } else {
        line <- line +
          # median line
          ggplot2::geom_line(
            stat = "summary",
            fun = "median",
            colour = mline.col,
            size = mline.size,
            ggplot2::aes(group = 1)
          )
      }
    } else {
      line <- line
    }

    # other themes
    line1 <- line +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::ylab("Normalized expression") + ggplot2::xlab("") +
      ggplot2::theme(
        axis.ticks.length = ggplot2::unit(0.1, "cm"),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          color = "black"
        ),
        strip.background = ggplot2::element_blank()
      ) +
      ggplot2::facet_wrap(~cluster_name, ncol = ncol, scales = "free")

    return(line1)
  } else {
    # ==========================================================================

    data <- data.frame(object$wide.res, check.names = FALSE)|>
      dplyr::arrange(as.numeric(as.character(cluster)))

    # prepare matrix
    if (object$type %in% c("mfuzz", "TCseq")) {
      mat <- data|>
        dplyr::arrange(as.numeric(as.character(cluster)))|>
        dplyr::select(-gene, -cluster, -membership)
    } else if (object$type == "wgcna") {
      mat <- data|>
        dplyr::arrange(as.numeric(as.character(cluster)))|>
        dplyr::select(-gene, -cluster, -modulecol)
    } else if (object$type == "scRNAdata") {
      mat <- data|>
        dplyr::arrange(as.numeric(as.character(cluster)))|>
        dplyr::select(-gene, -cluster)
    } else if (object$type == "monocle") {
      mat <- data|>
        dplyr::arrange(as.numeric(as.character(cluster)))|>
        dplyr::select(-gene, -cluster)
    } else {
      mat <- data|>
        dplyr::arrange(as.numeric(as.character(cluster)))|>
        dplyr::select(-gene, -cluster)
    }

    rownames(mat) <- data$gene

    if (object$geneMode == "all" | ncol(mat) > 20) {
      use_raster <- TRUE
    } else {
      use_raster <- FALSE
    }

    # sample orders
    if (!is.null(sample.order)) {
      mat <- mat[, sample.order]
    }

    # split info
    cl.info <- data.frame(table(data$cluster))|>
      dplyr::mutate(Var1 = as.numeric(as.character(Var1)))|>
      dplyr::arrange(Var1)
    cluster.num <- nrow(cl.info)

    subgroup <- lapply(seq_len(nrow(cl.info)), function(x) {
      nm <- rep(as.character(cl.info$Var1[x]), cl.info$Freq[x])
      paste("C", nm, sep = "")
    }) |> unlist()

    # cluster orders
    if (!is.null(cluster.order)) {
      subgroup <- factor(subgroup, levels = paste("C", cluster.order, sep = ""))
      cluster_row_slices <- FALSE
    } else {
      cluster_row_slices <- TRUE
    }

    # plot
    # =================== bar annotation for samples
    # sample group info
    if (object$geneMode == "all" & object$type == "scRNAdata") {
      # split info
      # celltype <- sapply(strsplit(colnames(mat),split = "\\|"), "[",2)
      celltype <- vapply(strsplit(colnames(mat), split = "\\|"), function(x) {
        x[2]
      }, character(1))
      cell.num.info <- table(celltype)[unique(celltype)]

      # order for column split
      if (is.null(sample.cell.order)) {
        column_split <- factor(rep(names(cell.num.info), cell.num.info),
                               levels = unique(celltype))
      } else {
        column_split <- factor(rep(names(cell.num.info), cell.num.info),
                               levels = sample.cell.order)
      }

      # assign colors for block
      if (is.null(sample.col)) {
        block.col <- seq_len(length(cell.num.info))
      } else {
        block.col <- sample.col
      }
    } else {
      if (is.null(sample.group)) {
        sample.info <- colnames(mat)

        # split columns
        if (is.null(HeatmapAnnotation)) {
          if (object$geneType == "branched") {
            if (ncol(mat) == 200) {
              column_split <- rep(c("branch1", "branch2"), each = 100)
            } else {
              column_split <- rep(levels(object$pseudotime),
                                  rep(100, ncol(mat) / 100))
            }
          } else {
            column_split <- NULL
          }
        } else {
          column_split <- column.split
        }
      } else {
        sample.info <- sample.group

        # split columns
        # column_split = sample.group
        column_split <- factor(sample.group, levels = unique(sample.group))
      }

      # assign colors for monocle input
      if (object$type != "monocle") {
        sample.info <- factor(sample.info, levels = unique(sample.info))

        # sample colors
        if (is.null(sample.col)) {
          # scol <- ggsci::pal_npg()(length(sample.info))

          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }

          scol <- circlize::rand_color(n = length(sample.info))
          names(scol) <- sample.info
        } else {
          scol <- sample.col
          names(scol) <- sample.info
        }
      } else {
        sample.info <- factor(object$pseudotime,
                              levels = unique(object$pseudotime))

        # sample colors
        if (is.null(pseudotime_col)) {
          if (object$geneType == "branched") {
            if (length(unique(object$pseudotime)) == 3) {
              pseudotime_col <- c("red", "grey80", "blue")
            } else {
              if (!requireNamespace("circlize", quietly = TRUE)) {
                stop("Package 'circlize' is required. Please install it.")
              }

              pseudotime_col <- circlize::rand_color(
                n = length(unique(object$pseudotime)))
            }
          } else {
            pseudotime_col <- c("blue", "red")
          }
        } else {
          pseudotime_col <- pseudotime_col
        }

        if (is.null(sample.col)) {
          if (object$type != "monocle") {
            if (!requireNamespace("circlize", quietly = TRUE)) {
              stop("Package 'circlize' is required. Please install it.")
            }

            scol <- circlize::rand_color(n = length(sample.info))
            names(scol) <- sample.info
          } else {
            if (object$geneType == "branched") {
              if (length(unique(object$pseudotime)) == 3) {
                scol <- rep(pseudotime_col,
                            table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              } else {
                scol <- rep(pseudotime_col,
                            table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
            } else {
              scol <- grDevices::colorRampPalette(pseudotime_col)(100)
              names(scol) <- sample.info
            }
          }
        } else {
          scol <- sample.col
          names(scol) <- sample.info
        }
      }
    }

    # top anno
    if (add.sampleanno == TRUE) {
      if (object$geneMode == "all" & object$type == "scRNAdata") {
        topanno <- ComplexHeatmap::HeatmapAnnotation(
          cluster = ComplexHeatmap::anno_block(
            gp = grid::gpar(fill = block.col),
            labels = NULL
          ),
          show_annotation_name = FALSE
        )
      } else {
        if (is.null(HeatmapAnnotation)) {
          topanno <- ComplexHeatmap::HeatmapAnnotation(
            sample = sample.info,
            col = list(sample = scol),
            gp = grid::gpar(col = ifelse(
              object$type == "monocle", NA, "white"
            )),
            show_legend = ifelse(object$type == "monocle", FALSE, TRUE),
            show_annotation_name = FALSE
          )
        } else {
          topanno <- HeatmapAnnotation
        }
      }
    } else {
      topanno <- NULL
    }

    # =================== bar annotation for clusters
    if (is.null(ctAnno.col)) {
      # if (requireNamespace("jjAnno", quietly = TRUE)){
      #   colanno <- jjAnno::useMyCol("stallion",n = cluster.num)
      # }else{
      #   stop("Package 'jjAnno' is required for this
      # functionality. Please install it.")
      # }
      if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required. Please install it.")
      }

      colanno <- circlize::rand_color(n = cluster.num)
    } else {
      colanno <- ctAnno.col
    }

    names(colanno) <- seq_len(cluster.num)
    # anno.block <- ComplexHeatmap::anno_block(
    # gp = grid::gpar(fill = colanno,col = NA),
    #                                          which = "row")

    align_to <- split(seq_len(nrow(mat)), subgroup)
    anno.block <- ComplexHeatmap::anno_block(
      align_to = align_to,
      panel_fun = function(index, nm) {
        npos <- as.numeric(unlist(strsplit(nm, split = "C"))[2])

        # rect
        grid::grid.rect(gp = grid::gpar(fill = colanno[npos], col = NA))

        # text
        if (annnoblock.text == TRUE) {
          grid::grid.text(
            label = paste("n:", length(index), sep = ""),
            rot = 90,
            gp = grid::gpar(
              col = annnoblock.gp[1],
              fontsize = as.numeric(annnoblock.gp[2])
            )
          )
        }
      },
      which = "row"
    )

    # =================== gene annotation for heatmap
    # whether mark your genes on plot
    if (!is.null(markGenes)) {
      # all genes
      rowGene <- rownames(mat)

      # tartget gene
      annoGene <- markGenes

      # add color for gene
      gene.col <- data|>
        dplyr::select(gene, cluster)|>
        dplyr::filter(gene %in% annoGene)

      purrr::map_df(seq_len(cluster.num), function(x) {
        tmp <- gene.col|>
          dplyr::filter(as.numeric(cluster) == x)|>
          dplyr::mutate(col = colanno[x])
      }) -> gene.col

      gene.col <- gene.col[match(annoGene, gene.col$gene), ]

      if (is.na(genes.gp[3])) {
        gcol <- gene.col$col
      } else {
        gcol <- genes.gp[3]
      }

      # get target gene index
      index <- match(annoGene, rowGene)

      # some genes annotation
      geneMark <- ComplexHeatmap::anno_mark(
        at = index,
        labels = annoGene,
        which = "row",
        side = markGenes.side,
        labels_gp = grid::gpar(
          fontface = genes.gp[1],
          fontsize = as.numeric(genes.gp[2]),
          col = gcol
        )
      )
    } else {
      geneMark <- NULL
    }

    # final annotation for heatmap
    right_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark,
                                                      cluster = anno.block)

    # =======================================================
    # return plot according to plot type
    if (object$type == "monocle" |
        object$geneMode == "all" | ncol(mat) > 20) {
      show_column_names <- FALSE
    } else {
      show_column_names <- TRUE
    }

    # legend for monocle heatmap
    if (object$geneType == "non-branched") {
      rg <- range(as.numeric(as.character(sample.info)))

      if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required. Please install it.")
      }

      col_fun2 <- circlize::colorRamp2(c(rg[1], rg[2]), pseudotime_col)
      lgd <- ComplexHeatmap::Legend(col_fun = col_fun2, title = "pseudotime")
      lgd_list <- list(lgd)
    } else if (object$geneType == "branched") {
      if (length(levels(sample.info)) == 3) {
        lgd <- ComplexHeatmap::Legend(
          labels = levels(sample.info),
          legend_gp = grid::gpar(fill = pseudotime_col),
          title = "branch"
        )
      } else {
        lgd <- ComplexHeatmap::Legend(
          labels = levels(sample.info),
          legend_gp = grid::gpar(fill = pseudotime_col),
          title = "branch"
        )
      }
      lgd_list <- list(lgd)
    } else {
      lgd_list <- NULL
    }

    # plot heatmap
    if (plot.type == "heatmap") {
      if (!is.null(row_annotation_obj)) {
        left_annotation_ht <- row_annotation_obj
      } else {
        left_annotation_ht <- NULL
      }

      # draw HT
      htf <-
        ComplexHeatmap::Heatmap(
          as.matrix(mat),
          name = "Z-score",
          cluster_columns = cluster_columns,
          show_row_names = show_row_names,
          border = border,
          column_split = column_split,
          row_split = subgroup,
          cluster_row_slices = cluster_row_slices,
          column_names_side = "top",
          show_column_names = show_column_names,
          # border = TRUE,
          top_annotation = topanno,
          left_annotation = left_annotation_ht,
          right_annotation = right_annotation,
          col = col_fun,
          use_raster = use_raster,
          ...
        )

      # draw
      ComplexHeatmap::draw(htf,
                           merge_legend = TRUE,
                           annotation_legend_list = lgd_list
      )
    } else {
      # ====================== heatmap + line
      rg <- range(mat)

      # ========================================================================
      # panel_fun for ggplot object
      # ========================================================================
      if (!is.null(gglist)) {
        anno_ggplot2 <- ComplexHeatmap::anno_zoom(
          align_to = align_to,
          which = "row",
          panel_fun = function(index, nm) {
            g <- gglist[[nm]]
            # g <- grid::grid.grabExpr(print(g))
            g <- grid::grid.grabExpr(grid::grid.draw(g))
            grid::pushViewport(grid::viewport())
            grid::grid.rect()
            grid::grid.draw(g)
            grid::popViewport()
          },
          size = grid::unit(as.numeric(ggplot.panel.arg[1]), "cm"),
          gap = grid::unit(as.numeric(ggplot.panel.arg[2]), "cm"),
          width = grid::unit(as.numeric(ggplot.panel.arg[3]), "cm"),
          side = "right",
          link_gp = grid::gpar(fill = ggplot.panel.arg[4],
                               col = ggplot.panel.arg[5])
        )
      } else {
        anno_ggplot2 <- NULL
      }

      # ========================================================================
      # panel_fun for line plot
      # ========================================================================
      panel_fun <- function(index, nm) {
        # whether add boxplot
        if (add.box == TRUE & add.line != TRUE) {
          xscale <- c(-0.1, 1.1)
        } else {
          # xscale = c(0,1)
          xscale <- c(-0.1, 1.1)
          panel_scale <- c(0.1, 0.9)
        }

        grid::pushViewport(grid::viewport(xscale = xscale, yscale = c(0, 1)))
        grid::grid.rect()

        # whether given multiple groups
        # if(is.null(mulGroup)){
        #   mulGroup <- ncol(mat)
        #
        #   # ================ calculate group columns index
        #   seqn <- data.frame(st = 1,sp = ncol(mat))
        # }else{
        #   mulGroup <- mulGroup
        #
        #   grid::grid.lines(x = c(0,1),y = rep(0.5,2),
        #                    gp = grid::gpar(col = "black",lty = "dashed"))
        #
        #   # ================ calculate group columns index
        #   cu <- cumsum(mulGroup)
        #   seqn <- data.frame(st = c(1,cu[1:(length(cu) - 1)] + 1),
        #                      sp = c(cu[1],cu[2:length(cu)]))
        # }

        if (object$geneMode == "all" & object$type == "scRNAdata") {
          mulGroup <- cell.num.info

          grid::grid.lines(
            x = c(0, 1),
            y = rep(0.5, 2),
            gp = grid::gpar(col = "black", lty = "dashed")
          )

          # ================ calculate group columns index
          cu <- cumsum(mulGroup)
          # seqn <- data.frame(st = c(1,cu[1:(length(cu) - 1)] + 1),
          #                    sp = c(cu[1],cu[2:length(cu)]))
          seqn <- data.frame(st = c(1, cu[seq(1, (length(cu) - 1), 1)] + 1),
                             sp = c(cu[1], cu[seq(2, length(cu), 1)]))
        } else {
          if (is.null(mulGroup)) {
            mulGroup <- ncol(mat)

            # ================ calculate group columns index
            seqn <- data.frame(st = 1, sp = ncol(mat))
          } else {
            mulGroup <- mulGroup

            grid::grid.lines(
              x = c(0, 1),
              y = rep(0.5, 2),
              gp = grid::gpar(col = "black", lty = "dashed")
            )

            # ================ calculate group columns index
            cu <- cumsum(mulGroup)
            # seqn <- data.frame(st = c(1,cu[1:(length(cu) - 1)] + 1),
            #                    sp = c(cu[1],cu[2:length(cu)]))
            seqn <- data.frame(st = c(1, cu[seq(1, length(cu) - 1, 1)] + 1),
                               sp = c(cu[1], cu[seq(2, length(cu), 1)]))
          }
        }


        # ======================================================================
        if (object$geneMode == "all" &&
            object$type == "scRNAdata") {
          # loop for multiple groups to create grobs
          purrr::map_dfr(seq_len(nrow(seqn)), function(x) {
            tmp <- seqn[x, ]

            # tmpmat <- mat[index, c(tmp$st:tmp$sp)]
            tmpmat <- mat[index, seq(tmp$st, tmp$sp, 1)]

            rg <- base::range(mat[index, ])

            # choose method
            if (set.md == "mean") {
              mdia <- base::mean(base::rowMeans(tmpmat))
            } else if (set.md == "median") {
              mdia <- stats::median(base::apply(tmpmat, 1, stats::median))
            } else {
              message("supply mean/median !")
            }

            res <- data.frame(x = x, val = mdia)
            return(res)
          }) -> cell.ave

          # lines grobs
          if (add.line == TRUE) {
            grid::grid.lines(
              x = scales::rescale(cell.ave$x, to = c(0, 1)),
              # y = scales::rescale(cell.ave$val,to = c(0,1),
              # from = c(rg[1] - 0.1,rg[2] + 0.1)),
              y = scales::rescale(cell.ave$val, to = c(0.1, 0.9)),
              gp = grid::gpar(lwd = 3, col = mline.col)
            )
          }
        } else {
          # ====================================================================
          # multiple lines
          base::lapply(seq_len(nrow(seqn)), function(x) {
            tmp <- seqn[x, ]
            # tmpmat <- mat[index, c(tmp$st:tmp$sp)]
            tmpmat <- mat[index, seq(tmp$st, tmp$sp, 1)]

            # choose method
            if (set.md == "mean") {
              mdia <- base::colMeans(tmpmat)
            } else if (set.md == "median") {
              mdia <- base::apply(tmpmat, 2, stats::median)
            } else {
              message("supply mean/median !")
            }

            # boxplot xpos
            pos <- scales::rescale(seq_len(ncol(tmpmat)), to = c(0, 1))

            # boxcol
            if (is.null(boxcol)) {
              boxcol <- rep("grey90", ncol(tmpmat))
            } else {
              boxcol <- boxcol
            }

            # boxplot grobs
            if (add.box == TRUE) {
              lapply(seq_len(ncol(tmpmat)), function(x) {
                ComplexHeatmap::grid.boxplot(
                  scales::rescale(
                    tmpmat[, x],
                    to = c(0, 1),
                    from = c(rg[1] - 0.5, rg[2] + 0.5)
                  ),
                  pos = pos[x],
                  direction = "vertical",
                  box_width = as.numeric(box.arg[1]),
                  outline = FALSE,
                  gp = grid::gpar(col = box.arg[2], fill = boxcol[x])
                )
              })
            }

            # points grobs
            if (add.point == TRUE) {
              grid::grid.points(
                x = scales::rescale(seq_len(ncol(
                  tmpmat
                )), to = panel_scale),
                y = scales::rescale(
                  mdia,
                  to = c(0, 1),
                  from = c(rg[1] - 0.5, rg[2] + 0.5)
                ),
                pch = as.numeric(point.arg[1]),
                gp = grid::gpar(fill = point.arg[2], col = point.arg[3]),
                size = grid::unit(as.numeric(point.arg[4]), "char")
              )
            }

            # lines grobs
            if (add.line == TRUE) {
              grid::grid.lines(
                x = scales::rescale(seq_len(ncol(
                  tmpmat
                )), to = panel_scale),
                y = scales::rescale(
                  mdia,
                  to = c(0, 1),
                  from = c(rg[1] - 0.5, rg[2] + 0.5)
                ),
                gp = grid::gpar(lwd = 3, col = mline.col[x])
              )
            }
          })
        }

        # get gene numbers
        grid.textbox <- utils::getFromNamespace("grid.textbox",
                                                "ComplexHeatmap")

        text <- paste("Gene Size:", nrow(mat[index, ]), sep = " ")
        grid.textbox(
          text,
          x = textbox.pos[1],
          y = textbox.pos[2],
          gp = grid::gpar(fontsize = textbox.size, fontface = "italic", ...)
        )

        grid::popViewport()
      }

      # whether annotate subgroups
      if (!is.null(subgroup.anno)) {
        align_to <- split(seq_len(nrow(mat)), subgroup)
        align_to <- align_to[subgroup.anno]
      } else {
        align_to <- subgroup
      }

      # anno link annotation
      anno <- ComplexHeatmap::anno_link(
        align_to = align_to,
        which = "row",
        panel_fun = panel_fun,
        size = grid::unit(as.numeric(panel.arg[1]), "cm"),
        gap = grid::unit(as.numeric(panel.arg[2]), "cm"),
        width = grid::unit(as.numeric(panel.arg[3]), "cm"),
        side = line.side,
        link_gp = grid::gpar(fill = panel.arg[4], col = panel.arg[5])
      )


      # ========================================================================
      # whether add go term annotations
      # ========================================================================
      if (!is.null(annoTerm.data)) {
        # load term info
        termanno <- annoTerm.data
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        } else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        } else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", "ratio")
        } else {
          message("No more than 4 columns!")
        }

        # term colors
        if (is.null(go.col)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }

          gocol <- circlize::rand_color(n = nrow(termanno))
        } else {
          gocol <- go.col
        }

        # term text size
        if (is.null(go.size)) {
          gosize <- rep(12, nrow(termanno))
        } else {
          if (go.size == "pval") {
            # loop for re-scaling pvalue
            purrr::map_df(unique(termanno$id), function(x) {
              tmp <- termanno|>
                dplyr::filter(id == x)|>
                dplyr::mutate(size = scales::rescale(-log10(pval),
                                                     to = term.text.limit))
            }) -> termanno.tmp

            gosize <- termanno.tmp$size
          } else {
            gosize <- go.size
          }
        }

        # add to termanno
        termanno <- termanno|>
          dplyr::ungroup()|>
          dplyr::mutate(col = gocol, fontsize = gosize)

        # to list
        lapply(seq_len(length(unique(
          termanno$id
        ))), function(x) {
          tmp <- termanno[which(termanno$id == unique(termanno$id)[x]), ]
          df <- data.frame(
            text = tmp$term,
            col = tmp$col,
            fontsize = tmp$fontsize
          )
          return(df)
        }) -> term.list

        # add names
        names(term.list) <- unique(termanno$id)

        # whether annotate subgroups
        if (!is.null(subgroup.anno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroup.anno]

          term.list <- term.list[subgroup.anno]
        } else {
          align_to2 <- subgroup
          term.list <- term.list
        }

        # textbox annotations
        # if(add.bar == TRUE){
        #   box.side = "left"
        # }else{
        #   box.side = "right"
        # }

        textbox <- ComplexHeatmap::anno_textbox(
          align_to2,
          term.list,
          word_wrap = word_wrap,
          add_new_line = add_new_line,
          side = annoTerm.mside,
          background_gp = grid::gpar(fill = termAnno.arg[1],
                                     col = termAnno.arg[2]),
          by = by.go
        )

        # final row annotation
        # if(line.side == "right"){
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(
        # cluster = anno.block,
        # line = anno,
        # textbox = textbox)
        #   left_annotation = NULL
        # }else{
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(
        # cluster = anno.block,
        # textbox = textbox)
        #   left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
        # }

        # GO bar anno function
        if (ncol(termanno) - 2 > 2) {
          anno_gobar <- function(data = NULL,
                                 bar.width = 0.1,
                                 # col = NA,
                                 align_to = NULL,
                                 panel.arg = panel.arg,
                                 ...) {
            # process data
            if (ncol(data) - 2 == 3) {
              data <- data|>
                dplyr::mutate(bary = -log10(pval))
            } else {
              data <- data|>
                dplyr::mutate(bary = ratio)
            }

            ComplexHeatmap::anno_zoom(
              align_to = align_to,
              which = "row",

              # =====================
              panel_fun = function(index, nm) {
                grid::pushViewport(grid::viewport(
                  xscale = c(0, 1),
                  yscale = c(0, 1)
                ))

                grid::grid.rect()

                # sub data
                tmp <- data|>
                  dplyr::filter(id == nm)
                # |>dplyr::arrange(bary)

                # bar grobs
                # grid::grid.rect(x = rep(0,nrow(tmp)),
                #                 y = scales::rescale(1:nrow(tmp),to = c(0,1)),
                #                 width = scales::rescale(tmp$log10P,
                # to = c(0,1)),
                #                 height = bar.width,
                #                 gp = grid::gpar(fill = tmp$col,col = col))

                grid::grid.segments(
                  x0 = rep(0, nrow(tmp)),
                  x1 = scales::rescale(rev(tmp$bary), to = c(0.1, 0.9)),
                  y0 = scales::rescale(seq_len(nrow(tmp)), to = c(0.1, 0.9)),
                  y1 = scales::rescale(seq_len(nrow(tmp)), to = c(0.1, 0.9)),
                  gp = grid::gpar(
                    lwd = bar.width,
                    col = rev(tmp$col),
                    lineend = "butt"
                  )
                )

                # add cluster name
                grid.textbox <- utils::getFromNamespace("grid.textbox",
                                                        "ComplexHeatmap")

                text <- nm
                grid.textbox(
                  text,
                  x = textbar.pos[1],
                  y = textbar.pos[2],
                  gp = grid::gpar(
                    fontsize = textbox.size,
                    fontface = "italic",
                    col = unique(tmp$col),
                    ...
                  )
                )

                grid::popViewport()
              },

              # =======================
              size = grid::unit(as.numeric(panel.arg[1]), "cm"),
              gap = grid::unit(as.numeric(panel.arg[2]), "cm"),
              width = grid::unit(as.numeric(panel.arg[3]), "cm"),
              side = "right",
              link_gp = grid::gpar(fill = termAnno.arg[1],
                                   col = termAnno.arg[2]),
              ...
            )
          }

          # ================================
          # bar anno
          baranno <- anno_gobar(
            data = termanno,
            align_to = align_to2,
            panel.arg = panel.arg,
            bar.width = bar.width
          )
        }
        # ======================================================================
        # whether add bar annotation
        # ======================================================================
        if (add.bar == TRUE) {
          baranno
        } else {
          baranno <- NULL
        }
      } else {
        # ======================================================
        # no GO annotation
        # if(line.side == "right"){
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(
        # cluster = anno.block,line = anno)
        #   left_annotation = NULL
        # }else{
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(
        # cluster = anno.block)
        #   left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
        # }
        textbox <- NULL
        baranno <- NULL
      }

      # ========================================================================
      # whether add kegg term annotations
      # ========================================================================
      if (!is.null(annoKegg.data)) {
        # load term info
        termanno <- annoKegg.data
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        } else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        } else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", "ratio")
        } else {
          message("No more than 4 columns!")
        }

        # term colors
        if (is.null(kegg.col)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }

          gocol <- circlize::rand_color(n = nrow(termanno))
        } else {
          gocol <- kegg.col
        }

        # term text size
        if (is.null(kegg.size)) {
          gosize <- rep(12, nrow(termanno))
        } else {
          if (kegg.size == "pval") {
            # loop for re-scaling pvalue
            purrr::map_df(unique(termanno$id), function(x) {
              tmp <- termanno|>
                dplyr::filter(id == x)|>
                dplyr::mutate(size = scales::rescale(-log10(pval),
                                                     to = term.text.limit))
            }) -> termanno.tmp

            gosize <- termanno.tmp$size
          } else {
            gosize <- kegg.size
          }
        }

        # add to termanno
        termanno <- termanno|>
          dplyr::ungroup()|>
          dplyr::mutate(col = gocol, fontsize = gosize)

        # to list
        lapply(seq_len(length(unique(
          termanno$id
        ))), function(x) {
          tmp <- termanno[which(termanno$id == unique(termanno$id)[x]), ]
          df <- data.frame(
            text = tmp$term,
            col = tmp$col,
            fontsize = tmp$fontsize
          )
          return(df)
        }) -> term.list

        # add names
        names(term.list) <- unique(termanno$id)

        # whether annotate subgroups
        if (!is.null(subgroup.anno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroup.anno]

          term.list <- term.list[subgroup.anno]
        } else {
          align_to2 <- subgroup
          term.list <- term.list
        }

        # anno_textbox
        textbox.kegg <- ComplexHeatmap::anno_textbox(
          align_to2,
          term.list,
          word_wrap = word_wrap,
          add_new_line = add_new_line,
          side = annoKegg.mside,
          background_gp = grid::gpar(fill = keggAnno.arg[1],
                                     col = keggAnno.arg[2]),
          by = by.kegg
        )

        # GO bar anno function
        if (ncol(termanno) - 2 > 2) {
          anno_keggbar <- function(data = NULL,
                                   bar.width = 0.1,
                                   # col = NA,
                                   align_to = NULL,
                                   panel.arg = panel.arg,
                                   ...) {
            # process data
            if (ncol(data) - 2 == 3) {
              data <- data|>
                dplyr::mutate(bary = -log10(pval))
            } else {
              data <- data|>
                dplyr::mutate(bary = ratio)
            }

            ComplexHeatmap::anno_zoom(
              align_to = align_to,
              which = "row",

              # =====================
              panel_fun = function(index, nm) {
                grid::pushViewport(grid::viewport(
                  xscale = c(0, 1),
                  yscale = c(0, 1)
                ))

                grid::grid.rect()

                # sub data
                tmp <- data|>
                  dplyr::filter(id == nm)
                # |>dplyr::arrange(bary)

                # bar grobs
                grid::grid.segments(
                  x0 = rep(0, nrow(tmp)),
                  x1 = scales::rescale(rev(tmp$bary), to = c(0.1, 0.9)),
                  y0 = scales::rescale(seq_len(nrow(tmp)), to = c(0.1, 0.9)),
                  y1 = scales::rescale(seq_len(nrow(tmp)), to = c(0.1, 0.9)),
                  gp = grid::gpar(
                    lwd = bar.width,
                    col = rev(tmp$col),
                    lineend = "butt"
                  )
                )

                # add cluster name
                grid.textbox <- utils::getFromNamespace("grid.textbox",
                                                        "ComplexHeatmap")

                text <- nm
                grid.textbox(
                  text,
                  x = textbar.pos[1],
                  y = textbar.pos[2],
                  gp = grid::gpar(
                    fontsize = textbox.size,
                    fontface = "italic",
                    col = unique(tmp$col),
                    ...
                  )
                )

                grid::popViewport()
              },

              # =======================
              size = grid::unit(as.numeric(panel.arg[1]), "cm"),
              gap = grid::unit(as.numeric(panel.arg[2]), "cm"),
              width = grid::unit(as.numeric(panel.arg[3]), "cm"),
              side = "right",
              link_gp = grid::gpar(fill = keggAnno.arg[1],
                                   col = keggAnno.arg[2]),
              ...
            )
          }

          # ================================
          # bar anno
          baranno.kegg <- anno_keggbar(
            data = termanno,
            align_to = align_to2,
            panel.arg = panel.arg,
            bar.width = bar.width
          )
        }
        # ======================================================================
        # whether add bar annotation
        # ======================================================================
        if (add.kegg.bar == TRUE) {
          baranno.kegg
        } else {
          baranno.kegg <- NULL
        }
      } else {
        textbox.kegg <- NULL
        baranno.kegg <- NULL
      }

      # ========================================================================
      # final row annotations
      # ========================================================================


      if (line.side == "right") {
        if (markGenes.side == "right") {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(
            gene = geneMark,
            cluster = anno.block,
            line = anno,
            anno_ggplot2 = anno_ggplot2,
            textbox = textbox,
            bar = baranno,
            textbox.kegg = textbox.kegg,
            baranno.kegg = baranno.kegg
          )

          if (!is.null(row_annotation_obj)) {
            left_annotation <- row_annotation_obj
          } else {
            left_annotation <- NULL
          }

          # left_annotation = NULL
        } else {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(
            cluster = anno.block,
            line = anno,
            anno_ggplot2 = anno_ggplot2,
            textbox = textbox,
            bar = baranno,
            textbox.kegg = textbox.kegg,
            baranno.kegg = baranno.kegg
          )

          if (!is.null(row_annotation_obj)) {
            left_annotation <- do.call(
              ComplexHeatmap::rowAnnotation,
              modifyList(list(gene = geneMark), row_annotation_obj)
            )
          } else {
            left_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark)
          }

          # left_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark)
        }
      } else {
        if (markGenes.side == "right") {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(
            gene = geneMark,
            cluster = anno.block,
            anno_ggplot2 = anno_ggplot2,
            textbox = textbox,
            bar = baranno,
            textbox.kegg = textbox.kegg,
            baranno.kegg = baranno.kegg
          )

          if (!is.null(row_annotation_obj)) {
            left_annotation <- do.call(
              ComplexHeatmap::rowAnnotation,
              modifyList(list(line = anno), row_annotation_obj)
            )
          } else {
            left_annotation <- ComplexHeatmap::rowAnnotation(line = anno)
          }

          # left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
        } else {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(
            cluster = anno.block,
            anno_ggplot2 = anno_ggplot2,
            textbox = textbox,
            bar = baranno,
            textbox.kegg = textbox.kegg,
            baranno.kegg = baranno.kegg
          )

          if (!is.null(row_annotation_obj)) {
            left_annotation <- do.call(
              ComplexHeatmap::rowAnnotation,
              modifyList(
                list(gene = geneMark, line = anno),
                row_annotation_obj
              )
            )
          } else {
            left_annotation <- ComplexHeatmap::rowAnnotation(line = anno,
                                                             gene = geneMark)
          }

          # left_annotation = ComplexHeatmap::rowAnnotation(line = anno,
          # gene = geneMark)
        }
      }

      # save
      if (object$type == "monocle" |
          object$geneMode == "all" | ncol(mat) > 20) {
        show_column_names <- FALSE
      } else {
        show_column_names <- TRUE
      }

      # pdf('test.pdf',height = 10,width = 10)
      htf <- ComplexHeatmap::Heatmap(
        as.matrix(mat),
        name = "Z-score",
        cluster_columns = cluster_columns,
        show_row_names = show_row_names,
        border = border,
        column_split = column_split,
        top_annotation = topanno,
        right_annotation = right_annotation2,
        left_annotation = left_annotation,
        column_names_side = "top",
        show_column_names = show_column_names,
        row_split = subgroup,
        cluster_row_slices = cluster_row_slices,
        col = col_fun,
        use_raster = use_raster,
        ...
      )

      # draw lines legend
      if (is.null(mulGroup)) {
        ComplexHeatmap::draw(htf,
                             merge_legend = TRUE,
                             annotation_legend_list = lgd_list
        )
      } else {
        if (is.null(lgd.label)) {
          lgd.label <- paste("group", seq_len(length(mulGroup)), sep = "")
        } else {
          lgd.label <- lgd.label
        }

        lgd_list2 <- ComplexHeatmap::Legend(
          labels = lgd.label,
          type = "lines",
          legend_gp = grid::gpar(col = mline.col, lty = 1)
        )

        if (!is.null(lgd_list)) {
          lgd_list_com <- ComplexHeatmap::packLegend(lgd_list, lgd_list2)
        } else {
          lgd_list_com <- lgd_list2
        }

        ComplexHeatmap::draw(htf,
                             annotation_legend_list = lgd_list_com,
                             merge_legend = TRUE
        )
      }
      # dev.off()
    }
  }
}
