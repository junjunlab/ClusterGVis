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
#' @param htColList list of heatmap col_range and col_color, default
#' list(col_range = c(-2, 0, 2),col_color = c("#08519C", "white", "#A50F15")).
#' @param border whether add border for heatmap, default TRUE.
#' @param plotType the plot type to choose which incuding "line","heatmap"
#' and "both".
#' @param msCol membership line color form Mfuzz cluster method results,
#' default c('#0099CC','grey90','#CC3333').
#' @param lineSize line size for line plot, default 0.1.
#' @param lineCol line color for line plot, default "grey90".
#' @param addMline whether add median line on plot, default TRUE.
#' @param mlineSize median line size, default 2.
#' @param mlineCol median line color, default "#CC3333".
#' @param ncol the columns for facet plot with line plot, default 4.
#' @param ctAnnoCol the heatmap cluster annotation bar colors, default NULL.
#' @param setMd the represent line method on heatmap-line plot(mean/median),
#'  default "median".
#' @param textboxPos the relative position of text in left-line plot,
#' default c(0.5,0.8).
#' @param textboxSize the text size of the text in left-line plot, default 8.
#' @param panelArg the settings for the left-line panel which are
#' panel size,gap,width,fill and col, default c(2,0.25,4,"grey90",NA).
#' @param ggplotPanelArg the settings for the ggplot2 object plot panel
#' which are
#' panel size,gap,width,fill and col, default c(2,0.25,4,"grey90",NA).
#' @param annoTermData the GO term annotation for the clusters, default NULL.
#' @param annoTermMside the wider GO term annotation box side, default "right".
#' @param termAnnoArg the settings for GO term panel annotations which
#' are fill and col,
#' default c("grey95","grey50").
#'
#' @param addBox whether add boxplot, default FALSE.
#' @param boxCol the box fill colors, default NULL.
#' @param boxArg this is related to boxplot width and border color, default
#' c(0.1,"grey50").
#' @param addPoint whether add point, default FALSE.
#' @param pointArg this is related to point shape,fill,color and size,
#' default c(19,"orange","orange",1).
#' @param addLine whether add line, default TRUE.
#' @param lineSide the line annotation side, default "right".
#'
#' @param markGenes the gene names to be added on plot, default NULL.
#' @param markGenesSide the gene label side, default "right".
#' @param genesGp gene labels graphics settings, default c('italic',10,NA).
#' @param goCol the GO term text colors, default NULL.
#' @param goSize the GO term text size(numeric or "pval"), default NULL.
#' @param mulGroup to draw multiple lines annotation, supply the groups numbers
#'  with vector, default NULL.
#' @param lgdLabel the lines annotation legend labels, default NULL.
#' @param showRowNames whether to show row names, default FALSE.
#' @param termTextLimit the GO term text size limit, default c(10,18).
#' @param subgroupAnno the sub-cluster for annotation, supply sub-cluster id,
#'  default NULL.
#' @param addBar whether add bar plot for GO enrichment, default FALSE.
#' @param barWidth the GO enrichment bar width, default 8.
#' @param textbarPos the barplot text relative position, default c(0.8,0.8).
#'
#' @param annnoblockText whether add cluster numbers on right block
#' annotation, default TRUE.
#' @param annnoblockGp right block annotation text color and size,
#' default c("white",8).
#' @param addSampleAnno whether add column annotation, default TRUE.
#' @param sampleGroup the column sample groups, default NULL.
#' @param sampleCol column annotation colors, default NULL.
#' @param sampleOrder the orders for column samples, default NULL.
#' @param heatmapAnnotation the 'heatmapAnnotation' object from 'ComplexHeatmap'
#' when you have multiple annotations, default NULL.
#' @param columnSplit how to split the columns when supply multiple column
#' annotations, default NULL.
#' @param clusterOrder the row cluster orders for user's own defination,
#' default NULL.
#' @param sampleCellOrder the celltype order when input is scRNA data and
#' "showAverage = FALSE"
#' for prepareDataFromscRNA.
#' @param annoKeggData the KEGG term annotation for the clusters, default NULL.
#' @param annoKeggMside the wider KEGG term annotation box side, default
#'  "right".
#' @param keggAnnoArg the settings for KEGG term panel annotations which are
#'  fill and col,
#' default c("grey95","grey50").
#' @param addKeggBar whether add bar plot for KEGG enrichment, default FALSE.
#' @param keggCol the KEGG term text colors, default NULL.
#' @param keggSize the KEGG term text size(numeric or "pval"), default NULL.
#' @param byGo the GO term text box style("anno_link" or "anno_block"),
#'  default "anno_link".
#' @param byKegg the KEGG term text box style("anno_link" or "anno_block"),
#'  default "anno_link".
#' @param wordWrap whether wrap the text, default TRUE.
#' @param addNewLine whether add new line when text is long, default TRUE.
#' @param clusterColumns whether cluster the columns, default FALSE.
#' @param pseudotimeCol the branch color control for monocle input data.
#' @param gglist a list of ggplot object to annotate each cluster, default NULL.
#' @param rowAnnotationObj Row annotation for heatmap, it is a
#'  `ComplexHeatmap::rowAnnotation()` object
#' when "markGenesSide" or ”lineSide“ is "right". Otherwise is a list of
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
#'   clusterMethod = "kmeans",
#'   clusterNum = 8
#' )
#'
#' # plot
#' visCluster(
#'   object = cm,
#'   plotType = "line"
#' )
#'
visCluster <- function(object = NULL,
                       # plot.data = NULL,
                       htColList = list(
                         col_range = c(-2, 0, 2),
                         col_color = c("#08519C", "white", "#A50F15")
                       ),
                       # ht.col = c("#08519C", "white", "#A50F15"),
                       border = TRUE,
                       plotType = c("line", "heatmap", "both"),
                       msCol = c("#0099CC", "grey90", "#CC3333"),
                       lineSize = 0.1,
                       lineCol = "grey90",
                       addMline = TRUE,
                       mlineSize = 2,
                       mlineCol = "#CC3333",
                       ncol = 4,
                       ctAnnoCol = NULL,
                       setMd = "median",
                       textboxPos = c(0.5, 0.8),
                       textboxSize = 8,
                       # panel size,gap,width,fill,col
                       panelArg = c(2, 0.25, 4, "grey90", NA),
                       ggplotPanelArg = c(2, 0.25, 4, "grey90", NA),
                       annoTermData = NULL,
                       annoTermMside = "right",
                       # textbox fill and col
                       termAnnoArg = c("grey95", "grey50"),
                       addBar = FALSE,
                       barWidth = 8,
                       textbarPos = c(0.8, 0.8),
                       goCol = NULL,
                       goSize = NULL,
                       byGo = "anno_link",
                       # KEGG term annotation
                       annoKeggData = NULL,
                       annoKeggMside = "right",
                       # textbox fill and col
                       keggAnnoArg = c("grey95", "grey50"),
                       addKeggBar = FALSE,
                       keggCol = NULL,
                       keggSize = NULL,
                       byKegg = "anno_link",
                       wordWrap = TRUE,
                       addNewLine = TRUE,
                       # boxplot,line.point annotation
                       addBox = FALSE,
                       boxCol = NULL,
                       # box with and border color
                       boxArg = c(0.1, "grey50"),
                       addPoint = FALSE,
                       # shape,fill,color,size
                       pointArg = c(19, "orange", "orange", 1),
                       addLine = TRUE,
                       lineSide = "right",
                       markGenes = NULL,
                       markGenesSide = "right",
                       genesGp = c("italic", 10, NA),
                       termTextLimit = c(10, 18),
                       mulGroup = NULL,
                       lgdLabel = NULL,
                       showRowNames = FALSE,
                       subgroupAnno = NULL,
                       annnoblockText = TRUE,
                       annnoblockGp = c("white", 8),
                       addSampleAnno = TRUE,
                       sampleGroup = NULL,
                       sampleCol = NULL,
                       sampleOrder = NULL,
                       clusterOrder = NULL,
                       sampleCellOrder = NULL,
                       heatmapAnnotation = NULL,
                       columnSplit = NULL,
                       clusterColumns = FALSE,
                       pseudotimeCol = NULL,
                       gglist = NULL,
                       rowAnnotationObj = NULL,
                       ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Please install it.")
  }

  ComplexHeatmap::ht_opt(message = FALSE)

  if (is.null(htColList[["col_range"]])) {
    col_range <- c(-2, 0, 2)
  } else {
    col_range <- htColList[["col_range"]]
  }

  if (is.null(htColList[["col_color"]])) {
    col_color <- c("#08519C", "white", "#A50F15")
  } else {
    col_color <- htColList[["col_color"]]
  }

  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Please install it.")
  }

  col_fun <- circlize::colorRamp2(col_range, col_color)

  plotType <- match.arg(plotType)

  # choose plot type
  if (plotType == "line") {
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
    if (!is.null(sampleOrder)) {
      data$cell_type <- factor(data$cell_type, levels = sampleOrder)
    }

    # basic plot
    line <-
      ggplot2::ggplot(data, ggplot2::aes(x = cell_type, y = norm_value))

    # type
    if (object$type %in% c("mfuzz", "TCseq")) {
      line <- line +
        ggplot2::geom_line(ggplot2::aes(color = membership, group = gene),
                           linewidth = lineSize) +
        ggplot2::scale_color_gradient2(
          low = msCol[1],
          mid = msCol[2],
          high = msCol[3],
          midpoint = 0.5
        )
    } else {
      line <- line +
        ggplot2::geom_line(ggplot2::aes(group = gene),
                           color = lineCol,
                           linewidth = lineSize
        )
    }

    if (addMline == TRUE) {
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
            linewidth = mlineSize,
            ggplot2::aes(group = 1, color = modulecol)
          ) +
          ggplot2::scale_color_manual(values = linec)
      } else {
        line <- line +
          # median line
          ggplot2::geom_line(
            stat = "summary",
            fun = "median",
            colour = mlineCol,
            linewidth = mlineSize,
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
    if (!is.null(sampleOrder)) {
      mat <- mat[, sampleOrder]
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
    if (!is.null(clusterOrder)) {
      subgroup <- factor(subgroup, levels = paste("C", clusterOrder, sep = ""))
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
      if (is.null(sampleCellOrder)) {
        column_split <- factor(rep(names(cell.num.info), cell.num.info),
                               levels = unique(celltype))
      } else {
        column_split <- factor(rep(names(cell.num.info), cell.num.info),
                               levels = sampleCellOrder)
      }

      # assign colors for block
      if (is.null(sampleCol)) {
        block.col <- seq_len(length(cell.num.info))
      } else {
        block.col <- sampleCol
      }
    } else {
      if (is.null(sampleGroup)) {
        sample.info <- colnames(mat)

        # split columns
        if (is.null(heatmapAnnotation)) {
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
          column_split <- columnSplit
        }
      } else {
        sample.info <- sampleGroup

        # split columns
        # column_split = sampleGroup
        column_split <- factor(sampleGroup, levels = unique(sampleGroup))
      }

      # assign colors for monocle input
      if (object$type != "monocle") {
        sample.info <- factor(sample.info, levels = unique(sample.info))

        # sample colors
        if (is.null(sampleCol)) {
          # scol <- ggsci::pal_npg()(length(sample.info))

          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }

          scol <- circlize::rand_color(n = length(sample.info))
          names(scol) <- sample.info
        } else {
          scol <- sampleCol
          names(scol) <- sample.info
        }
      } else {
        sample.info <- factor(object$pseudotime,
                              levels = unique(object$pseudotime))

        # sample colors
        if (is.null(pseudotimeCol)) {
          if (object$geneType == "branched") {
            if (length(unique(object$pseudotime)) == 3) {
              pseudotimeCol <- c("red", "grey80", "blue")
            } else {
              if (!requireNamespace("circlize", quietly = TRUE)) {
                stop("Package 'circlize' is required. Please install it.")
              }

              pseudotimeCol <- circlize::rand_color(
                n = length(unique(object$pseudotime)))
            }
          } else {
            pseudotimeCol <- c("blue", "red")
          }
        } else {
          pseudotimeCol <- pseudotimeCol
        }

        if (is.null(sampleCol)) {
          if (object$type != "monocle") {
            if (!requireNamespace("circlize", quietly = TRUE)) {
              stop("Package 'circlize' is required. Please install it.")
            }

            scol <- circlize::rand_color(n = length(sample.info))
            names(scol) <- sample.info
          } else {
            if (object$geneType == "branched") {
              if (length(unique(object$pseudotime)) == 3) {
                scol <- rep(pseudotimeCol,
                            table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              } else {
                scol <- rep(pseudotimeCol,
                            table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
            } else {
              scol <- grDevices::colorRampPalette(pseudotimeCol)(100)
              names(scol) <- sample.info
            }
          }
        } else {
          scol <- sampleCol
          names(scol) <- sample.info
        }
      }
    }

    # top anno
    if (addSampleAnno == TRUE) {
      if (object$geneMode == "all" & object$type == "scRNAdata") {
        topanno <- ComplexHeatmap::HeatmapAnnotation(
          cluster = ComplexHeatmap::anno_block(
            gp = grid::gpar(fill = block.col),
            labels = NULL
          ),
          show_annotation_name = FALSE
        )
      } else {
        if (is.null(heatmapAnnotation)) {
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
          topanno <- heatmapAnnotation
        }
      }
    } else {
      topanno <- NULL
    }

    # =================== bar annotation for clusters
    if (is.null(ctAnnoCol)) {
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
      colanno <- ctAnnoCol
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
        if (annnoblockText == TRUE) {
          grid::grid.text(
            label = paste("Num: ", length(index), sep = ""),
            rot = 90,
            gp = grid::gpar(
              col = annnoblockGp[1],
              fontsize = as.numeric(annnoblockGp[2])
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

      if (is.na(genesGp[3])) {
        gcol <- gene.col$col
      } else {
        gcol <- genesGp[3]
      }

      # get target gene index
      index <- match(annoGene, rowGene)

      # some genes annotation
      geneMark <- ComplexHeatmap::anno_mark(
        at = index,
        labels = annoGene,
        which = "row",
        side = markGenesSide,
        labels_gp = grid::gpar(
          fontface = genesGp[1],
          fontsize = as.numeric(genesGp[2]),
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

      col_fun2 <- circlize::colorRamp2(c(rg[1], rg[2]), pseudotimeCol)
      lgd <- ComplexHeatmap::Legend(col_fun = col_fun2, title = "pseudotime")
      lgd_list <- list(lgd)
    } else if (object$geneType == "branched") {
      if (length(levels(sample.info)) == 3) {
        lgd <- ComplexHeatmap::Legend(
          labels = levels(sample.info),
          legend_gp = grid::gpar(fill = pseudotimeCol),
          title = "branch"
        )
      } else {
        lgd <- ComplexHeatmap::Legend(
          labels = levels(sample.info),
          legend_gp = grid::gpar(fill = pseudotimeCol),
          title = "branch"
        )
      }
      lgd_list <- list(lgd)
    } else {
      lgd_list <- NULL
    }

    # plot heatmap
    if (plotType == "heatmap") {
      if (!is.null(rowAnnotationObj)) {
        left_annotation_ht <- rowAnnotationObj
      } else {
        left_annotation_ht <- NULL
      }

      # draw HT
      htf <-
        ComplexHeatmap::Heatmap(
          as.matrix(mat),
          name = "Z-score",
          cluster_columns = clusterColumns,
          show_row_names = showRowNames,
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
          size = grid::unit(as.numeric(ggplotPanelArg[1]), "cm"),
          gap = grid::unit(as.numeric(ggplotPanelArg[2]), "cm"),
          width = grid::unit(as.numeric(ggplotPanelArg[3]), "cm"),
          side = "right",
          link_gp = grid::gpar(fill = ggplotPanelArg[4],
                               col = ggplotPanelArg[5])
        )
      } else {
        anno_ggplot2 <- NULL
      }

      # ========================================================================
      # panel_fun for line plot
      # ========================================================================
      panel_fun <- function(index, nm) {
        # whether add boxplot
        if (addBox == TRUE & addLine != TRUE) {
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
            if (setMd == "mean") {
              mdia <- base::mean(base::rowMeans(tmpmat))
            } else if (setMd == "median") {
              mdia <- stats::median(base::apply(tmpmat, 1, stats::median))
            } else {
              message("supply mean/median !")
            }

            res <- data.frame(x = x, val = mdia)
            return(res)
          }) -> cell.ave

          # lines grobs
          if (addLine == TRUE) {
            grid::grid.lines(
              x = scales::rescale(cell.ave$x, to = c(0, 1)),
              # y = scales::rescale(cell.ave$val,to = c(0,1),
              # from = c(rg[1] - 0.1,rg[2] + 0.1)),
              y = scales::rescale(cell.ave$val, to = c(0.1, 0.9)),
              gp = grid::gpar(lwd = 3, col = mlineCol)
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
            if (setMd == "mean") {
              mdia <- base::colMeans(tmpmat)
            } else if (setMd == "median") {
              mdia <- base::apply(tmpmat, 2, stats::median)
            } else {
              message("supply mean/median !")
            }

            # boxplot xpos
            pos <- scales::rescale(seq_len(ncol(tmpmat)), to = c(0, 1))

            # boxCol
            if (is.null(boxCol)) {
              boxCol <- rep("grey90", ncol(tmpmat))
            } else {
              boxCol <- boxCol
            }

            # boxplot grobs
            if (addBox == TRUE) {
              lapply(seq_len(ncol(tmpmat)), function(x) {
                ComplexHeatmap::grid.boxplot(
                  scales::rescale(
                    tmpmat[, x],
                    to = c(0, 1),
                    from = c(rg[1] - 0.5, rg[2] + 0.5)
                  ),
                  pos = pos[x],
                  direction = "vertical",
                  box_width = as.numeric(boxArg[1]),
                  outline = FALSE,
                  gp = grid::gpar(col = boxArg[2], fill = boxCol[x])
                )
              })
            }

            # points grobs
            if (addPoint == TRUE) {
              grid::grid.points(
                x = scales::rescale(seq_len(ncol(
                  tmpmat
                )), to = c(0.1, 0.9)),
                y = scales::rescale(
                  mdia,
                  to = c(0, 1),
                  from = c(rg[1] - 0.5, rg[2] + 0.5)
                ),
                pch = as.numeric(pointArg[1]),
                gp = grid::gpar(fill = pointArg[2], col = pointArg[3]),
                size = grid::unit(as.numeric(pointArg[4]), "char")
              )
            }

            # lines grobs
            if (addLine == TRUE) {
              grid::grid.lines(
                x = scales::rescale(seq_len(ncol(
                  tmpmat
                )), to = c(0.1, 0.9)),
                y = scales::rescale(
                  mdia,
                  to = c(0, 1),
                  from = c(rg[1] - 0.5, rg[2] + 0.5)
                ),
                gp = grid::gpar(lwd = 3, col = mlineCol[x])
              )
            }
          })
        }

        # get gene numbers
        grid.textbox <- utils::getFromNamespace("grid.textbox",
                                                "ComplexHeatmap")

        text <- paste("Gene size:", nrow(mat[index, ]), sep = " ")
        grid.textbox(
          text,
          x = textboxPos[1],
          y = textboxPos[2],
          gp = grid::gpar(fontsize = textboxSize, fontface = "italic", ...)
        )

        grid::popViewport()
      }

      # whether annotate subgroups
      if (!is.null(subgroupAnno)) {
        align_to <- split(seq_len(nrow(mat)), subgroup)
        align_to <- align_to[subgroupAnno]
      } else {
        align_to <- subgroup
      }

      # anno link annotation
      anno <- ComplexHeatmap::anno_link(
        align_to = align_to,
        which = "row",
        panel_fun = panel_fun,
        size = grid::unit(as.numeric(panelArg[1]), "cm"),
        gap = grid::unit(as.numeric(panelArg[2]), "cm"),
        width = grid::unit(as.numeric(panelArg[3]), "cm"),
        side = lineSide,
        link_gp = grid::gpar(fill = panelArg[4], col = panelArg[5])
      )


      # ========================================================================
      # whether add go term annotations
      # ========================================================================
      if (!is.null(annoTermData)) {
        # load term info
        termanno <- annoTermData
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
        if (is.null(goCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }

          gocol <- circlize::rand_color(n = nrow(termanno))
        } else {
          gocol <- goCol
        }

        # term text size
        if (is.null(goSize)) {
          gosize <- rep(12, nrow(termanno))
        } else {
          if (goSize == "pval") {
            # loop for re-scaling pvalue
            purrr::map_df(unique(termanno$id), function(x) {
              tmp <- termanno|>
                dplyr::filter(id == x)|>
                dplyr::mutate(size = scales::rescale(-log10(pval),
                                                     to = termTextLimit))
            }) -> termanno.tmp

            gosize <- termanno.tmp$size
          } else {
            gosize <- goSize
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
        if (!is.null(subgroupAnno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroupAnno]

          term.list <- term.list[subgroupAnno]
        } else {
          align_to2 <- subgroup
          term.list <- term.list
        }

        # textbox annotations
        # if(addBar == TRUE){
        #   box.side = "left"
        # }else{
        #   box.side = "right"
        # }

        textbox <- ComplexHeatmap::anno_textbox(
          align_to2,
          term.list,
          word_wrap = wordWrap,
          add_new_line = addNewLine,
          side = annoTermMside,
          background_gp = grid::gpar(fill = termAnnoArg[1],
                                     col = termAnnoArg[2]),
          by = byGo
        )

        # final row annotation
        # if(lineSide == "right"){
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
                                 barWidth = 0.1,
                                 # col = NA,
                                 align_to = NULL,
                                 panelArg = panelArg,
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
                #                 height = barWidth,
                #                 gp = grid::gpar(fill = tmp$col,col = col))

                grid::grid.segments(
                  x0 = rep(0, nrow(tmp)),
                  x1 = scales::rescale(rev(tmp$bary), to = c(0.1, 0.9)),
                  y0 = scales::rescale(seq_len(nrow(tmp)), to = c(0.1, 0.9)),
                  y1 = scales::rescale(seq_len(nrow(tmp)), to = c(0.1, 0.9)),
                  gp = grid::gpar(
                    lwd = barWidth,
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
                  x = textbarPos[1],
                  y = textbarPos[2],
                  gp = grid::gpar(
                    fontsize = textboxSize,
                    fontface = "italic",
                    col = unique(tmp$col),
                    ...
                  )
                )

                grid::popViewport()
              },

              # =======================
              size = grid::unit(as.numeric(panelArg[1]), "cm"),
              gap = grid::unit(as.numeric(panelArg[2]), "cm"),
              width = grid::unit(as.numeric(panelArg[3]), "cm"),
              side = "right",
              link_gp = grid::gpar(fill = termAnnoArg[1],
                                   col = termAnnoArg[2]),
              ...
            )
          }

          # ================================
          # bar anno
          baranno <- anno_gobar(
            data = termanno,
            align_to = align_to2,
            panelArg = panelArg,
            barWidth = barWidth
          )
        }
        # ======================================================================
        # whether add bar annotation
        # ======================================================================
        if (addBar == TRUE) {
          baranno
        } else {
          baranno <- NULL
        }
      } else {
        # ======================================================
        # no GO annotation
        # if(lineSide == "right"){
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
      if (!is.null(annoKeggData)) {
        # load term info
        termanno <- annoKeggData
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
        if (is.null(keggCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }

          gocol <- circlize::rand_color(n = nrow(termanno))
        } else {
          gocol <- keggCol
        }

        # term text size
        if (is.null(keggSize)) {
          gosize <- rep(12, nrow(termanno))
        } else {
          if (keggSize == "pval") {
            # loop for re-scaling pvalue
            purrr::map_df(unique(termanno$id), function(x) {
              tmp <- termanno|>
                dplyr::filter(id == x)|>
                dplyr::mutate(size = scales::rescale(-log10(pval),
                                                     to = termTextLimit))
            }) -> termanno.tmp

            gosize <- termanno.tmp$size
          } else {
            gosize <- keggSize
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
        if (!is.null(subgroupAnno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroupAnno]

          term.list <- term.list[subgroupAnno]
        } else {
          align_to2 <- subgroup
          term.list <- term.list
        }

        # anno_textbox
        textbox.kegg <- ComplexHeatmap::anno_textbox(
          align_to2,
          term.list,
          word_wrap  = wordWrap,
          add_new_line  = addNewLine,
          side = annoKeggMside,
          background_gp = grid::gpar(fill = keggAnnoArg[1],
                                     col = keggAnnoArg[2]),
          by = byKegg
        )

        # GO bar anno function
        if (ncol(termanno) - 2 > 2) {
          anno_keggbar <- function(data = NULL,
                                   barWidth = 0.1,
                                   # col = NA,
                                   align_to = NULL,
                                   panelArg = panelArg,
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
                    lwd = barWidth,
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
                  x = textbarPos[1],
                  y = textbarPos[2],
                  gp = grid::gpar(
                    fontsize = textboxSize,
                    fontface = "italic",
                    col = unique(tmp$col),
                    ...
                  )
                )

                grid::popViewport()
              },

              # =======================
              size = grid::unit(as.numeric(panelArg[1]), "cm"),
              gap = grid::unit(as.numeric(panelArg[2]), "cm"),
              width = grid::unit(as.numeric(panelArg[3]), "cm"),
              side = "right",
              link_gp = grid::gpar(fill = keggAnnoArg[1],
                                   col = keggAnnoArg[2]),
              ...
            )
          }

          # ================================
          # bar anno
          baranno.kegg <- anno_keggbar(
            data = termanno,
            align_to = align_to2,
            panelArg = panelArg,
            barWidth = barWidth
          )
        }
        # ======================================================================
        # whether add bar annotation
        # ======================================================================
        if (addKeggBar == TRUE) {
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


      if (lineSide == "right") {
        if (markGenesSide == "right") {
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

          if (!is.null(rowAnnotationObj)) {
            left_annotation <- rowAnnotationObj
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

          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(
              ComplexHeatmap::rowAnnotation,
              modifyList(list(gene = geneMark), rowAnnotationObj)
            )
          } else {
            left_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark)
          }

          # left_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark)
        }
      } else {
        if (markGenesSide == "right") {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(
            gene = geneMark,
            cluster = anno.block,
            anno_ggplot2 = anno_ggplot2,
            textbox = textbox,
            bar = baranno,
            textbox.kegg = textbox.kegg,
            baranno.kegg = baranno.kegg
          )

          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(
              ComplexHeatmap::rowAnnotation,
              modifyList(list(line = anno), rowAnnotationObj)
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

          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(
              ComplexHeatmap::rowAnnotation,
              modifyList(
                list(gene = geneMark, line = anno),
                rowAnnotationObj
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
        cluster_columns = clusterColumns,
        show_row_names = showRowNames,
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
        if (is.null(lgdLabel)) {
          lgdLabel <- paste("group", seq_len(length(mulGroup)), sep = "")
        } else {
          lgdLabel <- lgdLabel
        }

        lgd_list2 <- ComplexHeatmap::Legend(
          labels = lgdLabel,
          type = "lines",
          legend_gp = grid::gpar(col = mlineCol, lty = 1)
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
