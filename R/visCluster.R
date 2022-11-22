#' @name visCluster
#' @author JunZhang
#' @title using visCluster to visualize cluster results from clusterData output
#'
#' @param object clusterData object, default NULL.
#' @param plot.type the plot type to choose which incuding "line","heatmap" and "both".
#' @param ms.col membership line color form Mfuzz cluster method results,
#' default c('#0099CC','grey90','#CC3333').
#' @param line.size line size for line plot, default 0.1.
#' @param line.col line color for line plot, default "grey90".
#' @param add.mline whether add median line on plot, default TRUE.
#' @param mline.size median line size, default 2.
#' @param mline.col median line color, default "#CC3333".
#' @param ncol the columns for facet plot with line plot, default 4.
#' @param ctAnno.col the heatmap cluster annotation bar colors, default NULL.
#' @param set.md the represent line method on heatmap-line plot(mean/median), default "median".
#' @param textbox.pos the relative position of text in left-line plot, default c(0.5,0.8).
#' @param textbox.size the text size of the text in left-line plot, default 8.
#' @param panel.arg the settings for the left-line panel which are
#' panel size,gap,width,fill and col, default c(2,0.25,4,"grey90",NA).
#' @param annoTerm.data the GO term annotation for the clusters, default NULL.
#' @param termAnno.arg the settings for GO term panel annotations which are fill and col,
#' default c("grey95","grey50").
#'
#' @param add.box whether add boxplot, default FALSE.
#' @param boxcol the box fill colors, default NULL.
#' @param box.arg this is related to boxplot width and border color, default c(0.1,"grey50").
#' @param add.point whether add point, default FALSE.
#' @param point.arg this is related to point shape,fill,color and size, default c(19,"orange","orange",1).
#' @param add.line whether add line, default TRUE.
#' @param line.side the line annotation side, default "right".
#'
#' @param markGenes the gene names to be added on plot, default NULL.
#' @param markGenes.side the gene label side, default "right".
#' @param genes.gp gene labels graphics settings, default c('italic',8,"black").
#' @param go.col the GO term text colors, default NULL.
#' @param go.size the GO term text size(numeric or "pval"), default NULL.
#' @param mulGroup to draw multipe lines annotation, supply the groups numbers with vector, default NULL.
#' @param lgd.label the lines annotation legend labels, default NULL.
#' @param show_row_names whether to show rownames, default FALSE.
#' @param term.text.limit the GO term text size limit, default c(5,18).
#'
#' @param ... othe aruguments passed by Heatmap fuction.
#'
#' @return a ggplot2 or Heatmap object.
#' @export
#'
#' @examples
#' \dontrun{
#' data("termanno")
#' data("exps")
#'
#' # mfuzz
#' cm <- clusterData(exp = exps,
#'                   cluster.method = "mfuzz",
#'                   cluster.num = 8)
#'
#' # plot
#' visCluster(object = cm,
#'            plot.type = "line")
#' }
globalVariables(c('cell_type', 'cluster.num', 'gene', 'membership', 'norm_value'))
visCluster <- function(object = NULL,
                       plot.type = c("line","heatmap","both"),
                       ms.col = c('#0099CC','grey90','#CC3333'),
                       line.size = 0.1,
                       line.col = "grey90",
                       add.mline = TRUE,
                       mline.size = 2,
                       mline.col = "#CC3333",
                       ncol = 4,
                       ctAnno.col = NULL,
                       set.md = "median",
                       textbox.pos = c(0.5,0.8),
                       textbox.size = 8,
                       # panel size,gap,width,fill,col
                       panel.arg = c(2,0.25,4,"grey90",NA),
                       annoTerm.data = NULL,
                       # textbox fill and col
                       termAnno.arg = c("grey95","grey50"),
                       add.box = FALSE,
                       boxcol = NULL,
                       # box with and border color
                       box.arg = c(0.1,"grey50"),
                       add.point = FALSE,
                       # shape,fill,color,size
                       point.arg = c(19,"orange","orange",1),
                       add.line = TRUE,
                       line.side = "right",
                       markGenes = NULL,
                       markGenes.side = "right",
                       genes.gp = c('italic',8,"black"),
                       go.col = NULL,
                       go.size = NULL,
                       term.text.limit = c(5,18),
                       mulGroup = NULL,
                       lgd.label = NULL,
                       show_row_names = FALSE,
                       ...){
  ComplexHeatmap::ht_opt(message = FALSE)
  plot.type <- match.arg(plot.type)

  # choose plot type
  if(plot.type == "line"){
    # process data
    data <- data.frame(object$long.res)

    # basic plot
    line <-
      ggplot2::ggplot(data,ggplot2::aes(x = cell_type,y = norm_value))

    # type
    if(object$type == "mfuzz"){
      line <- line +
        ggplot2::geom_line(ggplot2::aes(color = membership,group = gene),size = line.size) +
        ggplot2::scale_color_gradient2(low = ms.col[1],mid = ms.col[2],high = ms.col[3],
                                       midpoint = 0.5)

    }else{
      line <- line +
        ggplot2::geom_line(ggplot2::aes(group = gene),color = line.col,size = line.size)
    }

    if(add.mline == TRUE){
      line <- line +
        # median line
        ggplot2::geom_line(stat = "summary", fun = "median", colour = mline.col, size = mline.size,
                           ggplot2::aes(group = 1))
    }else{
      line <- line
    }

    # other themes
    line1 <- line +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::ylab('Normalized expression') + ggplot2::xlab('') +
      ggplot2::theme(axis.ticks.length = ggplot2::unit(0.1,'cm'),
                     axis.text.x = ggplot2::element_text(angle = 45,hjust = 1,color = 'black'),
                     strip.background = ggplot2::element_blank()) +
      ggplot2::facet_wrap(~cluster_name,ncol = ncol,scales = 'free')

    return(line1)

  }else{
    # ==========================================================================
    # process data
    data <- data.frame(object$wide.res)

    # prepare matrix
    if(object$type == "mfuzz"){
      mat <- data %>% dplyr::select(-gene,-cluster,-membership)
    }else{
      mat <- data %>% dplyr::select(-gene,-cluster)
    }

    rownames(mat) <- data$gene

    # split info
    cl.info <- data.frame(table(data$cluster))
    cluster.num <- nrow(cl.info)

    subgroup <- lapply(1:nrow(cl.info),function(x){
      nm <- rep(as.character(cl.info$Var1[x]),cl.info$Freq[x])
      paste("C",nm,sep = '')
    }) %>% unlist()

    # plot
    # =================== bar annotation for clusters
    if(is.null(ctAnno.col)){
      colanno <- jjAnno::useMyCol("stallion",n = cluster.num)
    }else{
      colanno <- ctAnno.col
    }

    names(colanno) <- 1:cluster.num
    anno.block <- ComplexHeatmap::anno_block(gp = grid::gpar(fill = colanno,col = NA),which = "row")

    # =================== gene annotation for heatmap
    # whether mark your genes on plot
    if(!is.null(markGenes)){
      # all genes
      rowGene <- rownames(mat)

      # tartget gene
      annoGene <- markGenes

      # get target gene index
      index <- match(annoGene,rowGene)

      # some genes annotation
      geneMark = gene = ComplexHeatmap::anno_mark(at = index,
                                                  labels = annoGene,
                                                  which = "row",
                                                  side = markGenes.side,
                                                  labels_gp = grid::gpar(fontface = genes.gp[1],
                                                                         fontsize = as.numeric(genes.gp[2]),
                                                                         col = genes.gp[3]))
    }else{
      geneMark = NULL
    }

    # final annotation for heatmap
    right_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark,cluster = anno.block)

    # =======================================================
    # return plot according to plot type
    if(plot.type == "heatmap"){
      # draw HT
      ComplexHeatmap::Heatmap(as.matrix(mat),
                              name = 'Z-score',
                              cluster_columns = FALSE,
                              show_row_names = FALSE,
                              row_split = subgroup,
                              column_names_side = "top",
                              border = TRUE,
                              right_annotation = right_annotation,
                              ...)
    }else{
      #====================== heatmap + line
      rg = range(mat)

      # # panel_fun for line plot
      # panel_fun = function(index, nm) {
      #   grid::pushViewport(grid::viewport(xscale = c(1,ncol(mat)), yscale = rg))
      #   grid::grid.rect()
      #
      #   # grid.xaxis(gp = gpar(fontsize = 8))
      #   # grid.annotation_axis(side = 'right',gp = gpar(fontsize = 8))
      #
      #   # choose method
      #   if(set.md == "mean"){
      #     mdia <- colMeans(mat[index, ])
      #   }else if(set.md == "median"){
      #     mdia <- apply(mat[index, ], 2, stats::median)
      #   }else{
      #     print("supply mean/median !")
      #   }
      #
      #   # get gene numbers
      #   text <- paste("Gene Size:",nrow(mat[index, ]),sep = ' ')
      #   ComplexHeatmap::grid.textbox(text,x = textbox.pos[1],y = textbox.pos[2],
      #                                gp = grid::gpar(fontsize = textbox.size,fontface = "italic"))
      #
      #   # grid.points(x = 1:ncol(m),y = mdia,
      #   #             pch = 19,
      #   #             gp = gpar(col = 'orange'))
      #
      #   grid::grid.lines(x = scales::rescale(1:ncol(mat),to = c(0,1)),
      #                    y = scales::rescale(mdia,to = c(0,1),from = rg),
      #                    gp = grid::gpar(lwd = 3,col = mline.col))
      #
      #   grid::popViewport()
      # }

      # ====================================================================
      # panel_fun for line plot
      panel_fun = function(index, nm) {

        # whether add boxplot
        if(add.box == TRUE & add.line != TRUE){
          xscale = c(-0.1,1.1)
        }else{
          xscale = c(0,1)
        }

        grid::pushViewport(grid::viewport(xscale = xscale, yscale = c(0,1)))
        grid::grid.rect()

        # grid.xaxis(gp = gpar(fontsize = 8))
        # grid.annotation_axis(side = 'right',gp = gpar(fontsize = 8))

        # # choose method
        # if(set.md == "mean"){
        #   mdia <- colMeans(mat[index, ])
        # }else if(set.md == "median"){
        #   mdia <- apply(mat[index, ], 2, stats::median)
        # }else{
        #   print("supply mean/median !")
        # }
        #
        # # boxplot xpos
        # pos = scales::rescale(1:ncol(mat),to = c(0,1))
        #
        # # boxcol
        # if(is.null(boxcol)){
        #   boxcol <- rep("grey90",ncol(mat))
        # }else{
        #   boxcol <- boxcol
        # }
        #
        # # boxplot grobs
        # if(add.box == TRUE){
        #   lapply(1:ncol(mat), function(x){
        #     ComplexHeatmap::grid.boxplot(scales::rescale(mat[index, ][,x],
        #                                                  to = c(0,1),
        #                                                  from = c(rg[1] - 0.5,rg[2] + 0.5)),
        #                                  pos = pos[x],
        #                                  direction = "vertical",
        #                                  box_width = as.numeric(box.arg[1]),
        #                                  outline = FALSE,
        #                                  gp = grid::gpar(col = box.arg[2],fill = boxcol[x]))
        #   })
        # }
        #
        # # points grobs
        # if(add.point == TRUE){
        #   grid::grid.points(x = scales::rescale(1:ncol(mat),to = c(0,1)),
        #                     y = scales::rescale(mdia,to = c(0,1),from = c(rg[1] - 0.5,rg[2] + 0.5)),
        #                     pch = as.numeric(point.arg[1]),
        #                     gp = grid::gpar(fill = point.arg[2],col = point.arg[3]),
        #                     size = grid::unit(as.numeric(point.arg[4]), "char"))
        # }
        #
        # # lines grobs
        # if(add.line == TRUE){
        #   grid::grid.lines(x = scales::rescale(1:ncol(mat),to = c(0,1)),
        #                    y = scales::rescale(mdia,to = c(0,1),from = c(rg[1] - 0.5,rg[2] + 0.5)),
        #                    gp = grid::gpar(lwd = 3,col = mline.col))
        # }

        # whether given multiple groups
        if(is.null(mulGroup)){
          mulGroup <- ncol(mat)

          # ================ calculate group columns index
          seqn <- data.frame(st = 1,
                             sp = ncol(mat))
        }else{
          mulGroup <- mulGroup

          grid::grid.lines(x = c(0,1),y = rep(0.5,2),
                           gp = grid::gpar(col = "black",lty = "dashed"))

          # ================ calculate group columns index
          cu <- cumsum(mulGroup)
          seqn <- data.frame(st = c(1,cu[1:(length(cu) - 1)] + 1),
                             sp = c(cu[1],cu[2:length(cu)]))
        }

        # loop for multiple groups to create grobs
        lapply(1:nrow(seqn), function(x){
          tmp <- seqn[x,]
          tmpmat <- mat[index, c(tmp$st:tmp$sp)]

          # choose method
          if(set.md == "mean"){
            mdia <- colMeans(tmpmat)
          }else if(set.md == "median"){
            mdia <- apply(tmpmat, 2, stats::median)
          }else{
            print("supply mean/median !")
          }

          # boxplot xpos
          pos = scales::rescale(1:ncol(tmpmat),to = c(0,1))

          # boxcol
          if(is.null(boxcol)){
            boxcol <- rep("grey90",ncol(tmpmat))
          }else{
            boxcol <- boxcol
          }

          # boxplot grobs
          if(add.box == TRUE){
            lapply(1:ncol(tmpmat), function(x){
              ComplexHeatmap::grid.boxplot(scales::rescale(tmpmat[,x],
                                                           to = c(0,1),
                                                           from = c(rg[1] - 0.5,rg[2] + 0.5)),
                                           pos = pos[x],
                                           direction = "vertical",
                                           box_width = as.numeric(box.arg[1]),
                                           outline = FALSE,
                                           gp = grid::gpar(col = box.arg[2],fill = boxcol[x]))
            })
          }

          # points grobs
          if(add.point == TRUE){
            grid::grid.points(x = scales::rescale(1:ncol(tmpmat),to = c(0,1)),
                              y = scales::rescale(mdia,to = c(0,1),from = c(rg[1] - 0.5,rg[2] + 0.5)),
                              pch = as.numeric(point.arg[1]),
                              gp = grid::gpar(fill = point.arg[2],col = point.arg[3]),
                              size = grid::unit(as.numeric(point.arg[4]), "char"))
          }

          # lines grobs
          if(add.line == TRUE){
            grid::grid.lines(x = scales::rescale(1:ncol(tmpmat),to = c(0,1)),
                             y = scales::rescale(mdia,to = c(0,1),from = c(rg[1] - 0.5,rg[2] + 0.5)),
                             gp = grid::gpar(lwd = 3,col = mline.col[x]))
          }
        })

        # get gene numbers
        grid.textbox <- utils::getFromNamespace("grid.textbox", "ComplexHeatmap")

        text <- paste("Gene Size:",nrow(mat[index, ]),sep = ' ')
        grid.textbox(text,x = textbox.pos[1],y = textbox.pos[2],
                     gp = grid::gpar(fontsize = textbox.size,
                                     fontface = "italic",
                                     ...))

        grid::popViewport()
      }

      # anno link annotation
      anno = ComplexHeatmap::anno_link(align_to = subgroup,
                                       which = "row",
                                       panel_fun = panel_fun,
                                       size = grid::unit(as.numeric(panel.arg[1]), "cm"),
                                       gap = grid::unit(as.numeric(panel.arg[2]), "cm"),
                                       width = grid::unit(as.numeric(panel.arg[3]), "cm"),
                                       side = line.side,
                                       link_gp = grid::gpar(fill = panel.arg[4],col = panel.arg[5]))

      # =====================================
      # whether add go term annotations
      if(!is.null(annoTerm.data)){
        # load term info
        termanno <- annoTerm.data
        if(ncol(termanno) == 2){
          colnames(termanno) <- c("id","term")
        }else if(ncol(termanno) == 3){
          colnames(termanno) <- c("id","term","pval")
        }else{
          print("No more than 3 columns!")
        }

        # term colors
        if(is.null(go.col)){
          gocol <- circlize::rand_color(n = nrow(termanno))
        }else{
          gocol <- go.col
        }

        # term text size
        if(is.null(go.size)){
          gosize <- rep(12,nrow(termanno))
        }else{
          if(go.size == "pval"){
            gosize <- scales::rescale(-log10(termanno$pval),to = term.text.limit)
          }else{
            gosize <- go.size
          }
        }

        # add to termanno
        termanno <- termanno %>%
          dplyr::mutate(col = gocol,fontsize = gosize)

        # to list
        lapply(1:length(unique(termanno$id)), function(x){
          tmp = termanno[which(termanno$id == unique(termanno$id)[x]),]
          df <- data.frame(text = tmp$term,
                           col = tmp$col,
                           fontsize = tmp$fontsize)
          return(df)
        }) -> term.list

        # add names
        names(term.list) <- unique(termanno$id)

        # textbox annotations
        textbox = ComplexHeatmap::anno_textbox(subgroup, term.list,
                                               word_wrap = TRUE,
                                               add_new_line = TRUE,
                                               background_gp = grid::gpar(fill = termAnno.arg[1],
                                                                          termAnno.arg[2]))

        # final row annotation
        # if(line.side == "right"){
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block,
        #                                                     line = anno,
        #                                                     textbox = textbox)
        #   left_annotation = NULL
        # }else{
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block,
        #                                                     textbox = textbox)
        #   left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
        # }

      }else{
        # ======================================================
        # no GO annotation
        # if(line.side == "right"){
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block,line = anno)
        #   left_annotation = NULL
        # }else{
        #   right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block)
        #   left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
        # }
        textbox = NULL
      }

      # ====================================================
      # final row annotations
      if(line.side == "right"){
        if(markGenes.side == "right"){
          right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark,
                                                            cluster = anno.block,
                                                            line = anno,
                                                            textbox = textbox)
          left_annotation = NULL
        }else{
          right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block,
                                                            line = anno,
                                                            textbox = textbox)
          left_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark)
        }

      }else{
        if(markGenes.side == "right"){
          right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark,
                                                            cluster = anno.block,
                                                            textbox = textbox)
          left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
        }else{
          right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block,
                                                            textbox = textbox)
          left_annotation = ComplexHeatmap::rowAnnotation(line = anno,
                                                          gene = geneMark)
        }
      }

      # save
      # pdf('test.pdf',height = 10,width = 10)
      htf <- ComplexHeatmap::Heatmap(as.matrix(mat),
                                     name = "Z-score",
                                     cluster_columns = FALSE,
                                     show_row_names = show_row_names,
                                     right_annotation = right_annotation2,
                                     left_annotation = left_annotation,
                                     column_names_side = "top",
                                     row_split = subgroup,
                                     ...)

      # draw lines legend
      if(is.null(mulGroup)){
        ComplexHeatmap::draw(htf)
      }else{
        if(is.null(lgd.label)){
          lgd.label <- paste("group",1:length(mulGroup),sep = '')
        }else{
          lgd.label <- lgd.label
        }

        lgd_list = list(
          ComplexHeatmap::Legend(labels = lgd.label,
                                 type = "lines",
                                 legend_gp = grid::gpar(col = mline.col, lty = 1)))

        ComplexHeatmap::draw(htf,annotation_legend_list = lgd_list,merge_legend = TRUE)
      }
      # dev.off()
    }

  }
}


###############################
#' This is a test data for this package
#' test data describtion
#'
#' @name termanno
#' @docType data
#' @author Junjun Lao
"termanno"

#' This is a test data for this package
#' test data describtion
#'
#' @name termanno2
#' @docType data
#' @author Junjun Lao
"termanno2"
