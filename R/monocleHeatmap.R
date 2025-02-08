globalVariables(c("Cluster", "Pseudotime", "State", "buildBranchCellDataSet",
                  "colorRampPalette", "fData", "genSmoothCurves", "pData", "vstExprs"))

#' traverseTree function
#'
#' @param g NULL
#' @param starting_cell NULL
#' @param end_cells NULL
#'
#'
traverseTree <- function(g, starting_cell, end_cells){

  if (requireNamespace("igraph", quietly = TRUE)) {
    # distance <- igraph::shortest.paths(g, v=starting_cell, to=end_cells)
    distance <- igraph::distances(g, v=starting_cell, to=end_cells)
    branchPoints <- which(igraph::degree(g) == 3)
    path <- igraph::shortest_paths(g, from = starting_cell, end_cells)
  }else{
    stop("Package 'igraph' is required for this functionality. Please install it.")
  }

  return(list(shortest_path = path$vpath, distance = distance,
              branch_points = intersect(branchPoints, unlist(path$vpath))))
}



#' Plots a pseudotime-ordered, row-centered heatmap which is slightly modified in monocle2
#'
#' @description The function plot_pseudotime_heatmap takes a CellDataSet object
#' (usually containing a only subset of significant genes) and generates smooth
#' expression curves much like plot_genes_in_pseudotime.
#' Then, it clusters these genes and plots them using the pheatmap package.
#' This allows you to visualize modules of genes that co-vary across pseudotime.
#'
#' @param cds_subset CellDataSet for the experiment (normally only the branching
#' genes detected with branchTest)
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' @param hclust_method The method used by pheatmap to perform hirearchical
#' clustering of the rows.
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param hmcols The color scheme for drawing the heatmap.
#' @param add_annotation_row Additional annotations to show for each row in the
#' heatmap. Must be a dataframe with one row for each row in the fData table of
#' cds_subset, with matching IDs.
#' @param add_annotation_col Additional annotations to show for each column in
#' the heatmap. Must be a dataframe with one row for each cell in the pData table
#' of cds_subset, with matching IDs.
#' @param show_rownames Whether to show the names for each row in the table.
#' @param use_gene_short_name Whether to use the short names for each row. If
#' FALSE, uses row IDs from the fData table.
#' @param scale_max The maximum value (in standard deviations) to show in the
#' heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the
#' heatmap. Values smaller than this are set to the min.
#' @param norm_method Determines how to transform expression values prior to
#' rendering
#' @param trend_formula A formula string specifying the model used in fitting
#' the spline curve for each gene/feature.
#' @param return_heatmap Whether to return the pheatmap object to the user.
#' @param cores Number of cores to use when smoothing the expression curves shown
#' in the heatmap.
#' @return A list of heatmap_matrix (expression matrix for the branch committment),
#' ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation
#' data.frame for the column).
#' @importFrom stats sd as.dist cor cutree
#' @export
plot_pseudotime_heatmap2 <- function(cds_subset,
                                     cluster_rows = TRUE,
                                     hclust_method = "ward.D2",
                                     num_clusters = 6,
                                     hmcols = NULL,
                                     add_annotation_row = NULL,
                                     add_annotation_col = NULL,
                                     show_rownames = FALSE,
                                     use_gene_short_name = TRUE,
                                     norm_method = c("log", "vstExprs"),
                                     scale_max = 3,
                                     scale_min = -3,
                                     trend_formula = '~sm.ns(Pseudotime, df=3)',
                                     return_heatmap = FALSE,
                                     cores = 1){
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(Biobase::pData(cds_subset)$Pseudotime), max(Biobase::pData(cds_subset)$Pseudotime),length.out = 100))

  if (requireNamespace("monocle", quietly = TRUE)) {
    m <- monocle::genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,
                                  relative_expr = T, new_data = newdata)


    #remove genes with no expression in any condition
    m=m[!apply(m,1,sum)==0,]

    norm_method <- match.arg(norm_method)

    # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
    if(norm_method == 'vstExprs' && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == FALSE) {
      m = monocle::vstExprs(cds_subset, expr_matrix=m)
    }
    else if(norm_method == 'log') {
      m = log10(m+pseudocount)
    }
  } else {
    warning("Cannot find monocle 'monocle' is not installed.")
  }



  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>scale_max] = scale_max
  m[m<scale_min] = scale_min

  heatmap_matrix <- m

  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1

  if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }else {
    bks <- seq(-3.1,3.1, length.out = length(hmcols))
  }

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ph <- pheatmap::pheatmap(heatmap_matrix,
                             useRaster = T,
                             cluster_cols=FALSE,
                             cluster_rows=cluster_rows,
                             show_rownames=F,
                             show_colnames=F,
                             clustering_distance_rows=row_dist,
                             clustering_method = hclust_method,
                             cutree_rows=num_clusters,
                             silent=TRUE,
                             filename=NA,
                             breaks=bks,
                             border_color = NA,
                             color=hmcols)
  } else {
    warning("Cannot create heatmap. 'pheatmap' is not installed.")
  }



  if(cluster_rows) {
    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  } else {
    annotation_row <- NULL
  }

  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(Biobase::fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }

  if(!is.null(add_annotation_col)) {
    if(nrow(add_annotation_col) != 100) {
      stop('add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!')
    }
    annotation_col <- add_annotation_col
  } else {
    annotation_col <- NA
  }

  if (use_gene_short_name == TRUE) {
    if (is.null(Biobase::fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(Biobase::fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)

      row_ann_labels <- as.character(Biobase::fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }else {
    feature_label <- row.names(heatmap_matrix)
    if(!is.null(annotation_row))
      row_ann_labels <- row.names(annotation_row)
  }

  row.names(heatmap_matrix) <- feature_label
  if(!is.null(annotation_row))
    row.names(annotation_row) <- row_ann_labels

  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ph_res <- pheatmap::pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                                 useRaster = T,
                                 cluster_cols = FALSE,
                                 cluster_rows = cluster_rows,
                                 show_rownames=show_rownames,
                                 show_colnames=F,
                                 #scale="row",
                                 clustering_distance_rows=row_dist, #row_dist
                                 clustering_method = hclust_method, #ward.D2
                                 cutree_rows=num_clusters,
                                 # cutree_cols = 2,
                                 annotation_row=annotation_row,
                                 annotation_col=annotation_col,
                                 treeheight_row = 20,
                                 breaks=bks,
                                 fontsize = 6,
                                 color=hmcols,
                                 border_color = NA,
                                 silent=TRUE,
                                 filename=NA
    )

  } else {
    warning("Cannot create heatmap. 'pheatmap' is not installed.")
  }


  # ============================================================================
  # prepare data
  wide.res <- cbind(heatmap_matrix,annotation_row) %>%
    data.frame(.,check.names = FALSE) %>%
    dplyr::mutate(gene = rownames(.),.before = 1) %>%
    dplyr::rename(cluster = Cluster)

  # wide to long
  df <- reshape2::melt(wide.res,
                       id.vars = c('cluster','gene'),
                       variable.name = 'cell_type',
                       value.name = 'norm_value') %>%
    dplyr::mutate(cell_type = as.numeric(as.character(cell_type)))

  # add cluster name
  df$cluster_name <- paste('cluster ',df$cluster,sep = '')

  # add gene number
  cl.info <- data.frame(table(wide.res$cluster)) %>%
    dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>%
    dplyr::arrange(Var1)

  id <- unique(df$cluster_name)
  purrr::map_df(seq_along(id),function(x){
    tmp <- df %>%
      dplyr::filter(cluster_name == id[x])

    tmp %>%
      dplyr::mutate(cluster_name = paste(cluster_name," (",cl.info$Freq[x],")",sep = ''))
  }) -> df

  # cluster order
  df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",cl.info$Var1,
                                                           " (",cl.info$Freq,")",sep = ''))

  # return
  prepared_data <- list(wide.res = wide.res,
                        long.res = df,
                        type = "monocle",
                        geneMode = "all",
                        geneType = "non-branched",
                        pseudotime = newdata$Pseudotime)

  if (return_heatmap == TRUE){
    grid::grid.rect(gp=grid::gpar("fill", col=NA))
    grid::grid.draw(ph_res$gtable)
    return(ph_res)
  }else{
    return(prepared_data)
  }
}


#'  Create a heatmap to demonstrate the bifurcation of gene expression along two
#'  branchs which is slightly modified in monocle2
#'
#'  @description returns a heatmap that shows changes in both lineages at the same
#'  time.
#'  It also requires that you choose a branch point to inspect.
#'  Columns are points in pseudotime, rows are genes, and the beginning of
#'  pseudotime is in the middle of the heatmap.
#'  As you read from the middle of the heatmap to the right, you are following
#'  one lineage through pseudotime. As you read left, the other.
#'  The genes are clustered hierarchically, so you can visualize modules of genes
#'  that have similar lineage-dependent expression patterns.
#'
#' @param cds_subset CellDataSet for the experiment (normally only the branching
#' genes detected with branchTest)
#' @param branch_point The ID of the branch point to visualize. Can only be used
#' when reduceDimension is called with method = "DDRTree".
#' @param branch_states The two states to compare in the heatmap. Mutually
#' exclusive with branch_point.
#' @param branch_labels The labels for the branchs.
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' @param hclust_method The method used by pheatmap to perform hirearchical
#' clustering of the rows.
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param hmcols The color scheme for drawing the heatmap.
#' @param branch_colors The colors used in the annotation strip indicating the
#' pre- and post-branch cells.
#' @param add_annotation_row Additional annotations to show for each row in the
#' heatmap. Must be a dataframe with one row for each row in the fData table of
#' cds_subset, with matching IDs.
#' @param add_annotation_col Additional annotations to show for each column in
#' the heatmap. Must be a dataframe with one row for each cell in the pData table
#' of cds_subset, with matching IDs.
#' @param show_rownames Whether to show the names for each row in the table.
#' @param use_gene_short_name Whether to use the short names for each row. If
#' FALSE, uses row IDs from the fData table.
#' @param scale_max The maximum value (in standard deviations) to show in the
#' heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the
#' heatmap. Values smaller than this are set to the min.
#' @param norm_method Determines how to transform expression values prior to
#' rendering
#' @param trend_formula A formula string specifying the model used in fitting
#' the spline curve for each gene/feature.
#' @param return_heatmap Whether to return the pheatmap object to the user.
#' @param cores Number of cores to use when smoothing the expression curves shown
#' in the heatmap.
#' @param ... Additional arguments passed to buildBranchCellDataSet
#' @return A list of heatmap_matrix (expression matrix for the branch committment),
#' ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation
#' data.frame for the column).
#' @importFrom stats sd as.dist cor cutree
#' @export
plot_genes_branched_heatmap2 <- function(cds_subset = NULL,
                                         branch_point = 1,
                                         branch_states = NULL,
                                         branch_labels = c("Cell fate 1", "Cell fate 2"),
                                         cluster_rows = TRUE,
                                         hclust_method = "ward.D2",
                                         num_clusters = 6,
                                         hmcols = NULL,
                                         branch_colors = c('#979797', '#F05662', '#7990C8'),
                                         add_annotation_row = NULL,
                                         add_annotation_col = NULL,
                                         show_rownames = FALSE,
                                         use_gene_short_name = TRUE,
                                         scale_max = 3,
                                         scale_min = -3,
                                         norm_method = c("log", "vstExprs"),
                                         trend_formula = '~sm.ns(Pseudotime, df=3) * Branch',
                                         return_heatmap = FALSE,
                                         cores = 1, ...) {
  if (requireNamespace("monocle", quietly = TRUE)) {
    cds <- NA
    new_cds <- monocle::buildBranchCellDataSet(cds_subset,
                                               branch_states=branch_states,
                                               branch_point=branch_point,
                                               progenitor_method = 'duplicate',
                                               ...)
  } else {
    warning("Cannot find monocle 'monocle' is not installed.")
  }

  new_cds@dispFitInfo <- cds_subset@dispFitInfo

  if(is.null(branch_states)) {
    progenitor_state <- subset(Biobase::pData(cds_subset), Pseudotime == 0)[, 'State']
    branch_states <- setdiff(Biobase::pData(cds_subset)$State, progenitor_state)
  }

  col_gap_ind <- 101
  # newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
  # newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100))

  newdataA <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), Branch = as.factor(unique(as.character(Biobase::pData(new_cds)$Branch))[1]))
  newdataB <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), Branch = as.factor(unique(as.character(Biobase::pData(new_cds)$Branch))[2]))

  if (requireNamespace("monocle", quietly = TRUE)) {
    BranchAB_exprs <- monocle::genSmoothCurves(new_cds[, ], cores=cores, trend_formula = trend_formula,
                                               relative_expr = T, new_data = rbind(newdataA, newdataB))
  } else {
    warning("Cannot find monocle 'monocle' is not installed.")
  }


  BranchA_exprs <- BranchAB_exprs[, 1:100]
  BranchB_exprs <- BranchAB_exprs[, 101:200]

  #common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
  common_ancestor_cells <- row.names(Biobase::pData(new_cds)[Biobase::pData(new_cds)$State == setdiff(Biobase::pData(new_cds)$State, branch_states),])
  BranchP_num <- (100 - floor(max(Biobase::pData(new_cds)[common_ancestor_cells, 'Pseudotime'])))
  BranchA_num <- floor(max(Biobase::pData(new_cds)[common_ancestor_cells, 'Pseudotime']))
  BranchB_num <- BranchA_num

  norm_method <- match.arg(norm_method)

  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if (requireNamespace("monocle", quietly = TRUE)) {
    if(norm_method == 'vstExprs') {
      BranchA_exprs <- monocle::vstExprs(new_cds, expr_matrix=BranchA_exprs)
      BranchB_exprs <- monocle::vstExprs(new_cds, expr_matrix=BranchB_exprs)
    }else if(norm_method == 'log') {
      BranchA_exprs <- log10(BranchA_exprs + 1)
      BranchB_exprs <- log10(BranchB_exprs + 1)
    }
  } else {
    warning("Cannot find monocle 'monocle' is not installed.")
  }


  heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], BranchB_exprs)

  heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
  heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
  heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE,]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix>scale_max] = scale_max
  heatmap_matrix[heatmap_matrix<scale_min] = scale_min

  heatmap_matrix_ori <- heatmap_matrix
  heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch

  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1

  exp_rng <- range(heatmap_matrix) #bks is based on the expression range
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
  if(is.null(hmcols)) {
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }

  # prin  t(hmcols)
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ph <- pheatmap::pheatmap(heatmap_matrix,
                             useRaster = T,
                             cluster_cols=FALSE,
                             cluster_rows=TRUE,
                             show_rownames=F,
                             show_colnames=F,
                             #scale="row",
                             clustering_distance_rows=row_dist,
                             clustering_method = hclust_method,
                             cutree_rows=num_clusters,
                             silent=TRUE,
                             filename=NA,
                             breaks=bks,
                             color=hmcols
                             #color=hmcols#,
                             # filename="expression_pseudotime_pheatmap.pdf",
    )
  } else {
    warning("Cannot create heatmap. 'pheatmap' is not installed.")
  }


  #save(heatmap_matrix, row_dist, num_clusters, hmcols, ph, branchTest_df, qval_lowest_thrsd, branch_labels, BranchA_num, BranchP_num, BranchB_num, file = 'heatmap_matrix')

  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))

  if(!is.null(add_annotation_row)) {
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
    # annotation_row$bif_time <- add_annotation_row[as.character(Biobase::fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }

  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), "Cell Type" = c(rep(branch_labels[1], BranchA_num),
                                                                                      rep("Pre-branch",  2 * BranchP_num),
                                                                                      rep(branch_labels[2], BranchB_num)))

  colnames(annotation_col) <- "Cell Type"

  if(!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, add_annotation_col[Biobase::fData(cds[row.names(annotation_col), ])$gene_short_name, 1])
  }

  names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])

  annotation_colors=list("Cell Type"=branch_colors)

  names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)

  if (use_gene_short_name == TRUE) {
    if (is.null(Biobase::fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(Biobase::fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)

      row_ann_labels <- as.character(Biobase::fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }

  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ph_res <- pheatmap::pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                                 useRaster = T,
                                 cluster_cols=FALSE,
                                 cluster_rows=TRUE,
                                 show_rownames=show_rownames,
                                 show_colnames=F,
                                 #scale="row",
                                 clustering_distance_rows=row_dist, #row_dist
                                 clustering_method = hclust_method, #ward.D2
                                 cutree_rows=num_clusters,
                                 # cutree_cols = 2,
                                 annotation_row=annotation_row,
                                 annotation_col=annotation_col,
                                 annotation_colors=annotation_colors,
                                 gaps_col = col_gap_ind,
                                 treeheight_row = 20,
                                 breaks=bks,
                                 fontsize = 6,
                                 color=hmcols,
                                 border_color = NA,
                                 silent=TRUE)
  } else {
    warning("Cannot create heatmap. 'pheatmap' is not installed.")
  }



  # ============================================================================
  # prepare data
  wide.res <- cbind(heatmap_matrix,annotation_row) %>%
    data.frame(.,check.names = FALSE) %>%
    dplyr::mutate(gene = rownames(.),.before = 1) %>%
    dplyr::rename(cluster = Cluster)

  # wide to long
  df <- reshape2::melt(wide.res,
                       id.vars = c('cluster','gene'),
                       variable.name = 'cell_type',
                       value.name = 'norm_value') %>%
    dplyr::mutate(cell_type = as.numeric(as.character(cell_type)))

  # add cluster name
  df$cluster_name <- paste('cluster ',df$cluster,sep = '')

  # add gene number
  cl.info <- data.frame(table(wide.res$cluster)) %>%
    dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>%
    dplyr::arrange(Var1)

  id <- unique(df$cluster_name)
  purrr::map_df(seq_along(id),function(x){
    tmp <- df %>%
      dplyr::filter(cluster_name == id[x])

    tmp %>%
      dplyr::mutate(cluster_name = paste(cluster_name," (",cl.info$Freq[x],")",sep = ''))
  }) -> df

  # cluster order
  df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",cl.info$Var1,
                                                           " (",cl.info$Freq,")",sep = ''))

  # return
  prepared_data <- list(wide.res = wide.res,
                        long.res = df,
                        type = "monocle",
                        geneMode = "all",
                        geneType = "branched",
                        pseudotime = annotation_col$`Cell Type`)

  if (return_heatmap == TRUE){
    grid::grid.rect(gp=grid::gpar("fill", col=NA))
    grid::grid.draw(ph_res$gtable)
    return(ph_res)
    # return(list(BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs, heatmap_matrix = heatmap_matrix,
    #             heatmap_matrix_ori = heatmap_matrix_ori, ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, hmcols = hmcols,
    #             annotation_colors = annotation_colors, annotation_row = annotation_row, annotation_col = annotation_col,
    #             ph_res = ph_res))
  }else{
    return(prepared_data)
  }
}


#  Modified function: Plot heatmap of 3 branches with the same coloring. Each CDS
#  subset has to have the same set of genes.
#' Create a heatmap to demonstrate the bifurcation of gene expression along
#' multiple branches
#'
#' @param cds CellDataSet for the experiment (normally only the branching genes
#' detected with BEAM)
#' @param branches The terminal branches (states) on the developmental tree you
#' want to investigate.
#' @param branches_name Name (for example, cell type) of branches you believe the
#' cells on the branches are associated with.
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' @param hclust_method The method used by pheatmap to perform hirearchical
#' clustering of the rows.
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param hmcols The color scheme for drawing the heatmap.
#' @param add_annotation_row Additional annotations to show for each row in the
#' heatmap. Must be a dataframe with one row for each row in the fData table of
#' cds_subset, with matching IDs.
#' @param add_annotation_col Additional annotations to show for each column in
#' the heatmap. Must be a dataframe with one row for each cell in the pData
#' table of cds_subset, with matching IDs.
#' @param show_rownames Whether to show the names for each row in the table.
#' @param use_gene_short_name Whether to use the short names for each row. If
#' FALSE, uses row IDs from the fData table.
#' @param norm_method Determines how to transform expression values prior to
#' rendering
#' @param scale_max The maximum value (in standard deviations) to show in the
#' heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the
#' heatmap. Values smaller than this are set to the min.
#' @param trend_formula A formula string specifying the model used in fitting
#' the spline curve for each gene/feature.
#' @param return_heatmap Whether to return the pheatmap object to the user.
#' @param cores Number of cores to use when smoothing the expression curves shown
#' in the heatmap.
#' @return A list of heatmap_matrix (expression matrix for the branch committment),
#' ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation
#' data.frame for the column).
#' @export
plot_multiple_branches_heatmap2 <- function(cds = NULL,
                                            branches,
                                            branches_name = NULL,
                                            cluster_rows = TRUE,
                                            hclust_method = "ward.D2",
                                            num_clusters = 6,
                                            hmcols = NULL,
                                            add_annotation_row = NULL,
                                            add_annotation_col = NULL,
                                            show_rownames = FALSE,
                                            use_gene_short_name = TRUE,
                                            norm_method = c("vstExprs", "log"),
                                            scale_max=3,
                                            scale_min=-3,
                                            trend_formula = '~sm.ns(Pseudotime, df=3)',
                                            return_heatmap=FALSE,
                                            cores=1){
  pseudocount <- 1
  if(!(all(branches %in% Biobase::pData(cds)$State)) & length(branches) == 1){
    stop('This function only allows to make multiple branch plots where branches is included in the pData')
  }

  branch_label <- branches
  if(!is.null(branches_name)){
    if(length(branches) != length(branches_name)){
      stop('branches_name should have the same length as branches')
    }
    branch_label <- branches_name
  }

  #test whether or not the states passed to branches are true branches (not truncks) or there are terminal cells
  g <- cds@minSpanningTree
  m <- NULL
  # branche_cell_num <- c()
  for(branch_in in branches) {
    branches_cells <- row.names(subset(Biobase::pData(cds), State == branch_in))
    root_state <- subset(Biobase::pData(cds), Pseudotime == 0)[, 'State']
    root_state_cells <- row.names(subset(Biobase::pData(cds), State == root_state))

    if(cds@dim_reduce_type != 'ICA') {
      root_state_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_state_cells, ], sep = ''))
      branches_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[branches_cells, ], sep = ''))
    }
    root_cell <- root_state_cells[which(igraph::degree(g, v = root_state_cells) == 1)]
    tip_cell <- branches_cells[which(igraph::degree(g, v = branches_cells) == 1)]

    traverse_res <- traverseTree(g, root_cell, tip_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])

    if(cds@dim_reduce_type != 'ICA') {
      pc_ind <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
      path_cells <- row.names(pc_ind)[paste('Y_', pc_ind[, 1], sep = '') %in% path_cells]
    }

    cds_subset <- cds[, path_cells]

    newdata <- data.frame(Pseudotime = seq(0, max(Biobase::pData(cds_subset)$Pseudotime),length.out = 100))

    if (requireNamespace("monocle", quietly = TRUE)) {
      tmp <- monocle::genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,
                                      relative_expr = T, new_data = newdata)
    } else {
      warning("Cannot find monocle 'monocle' is not installed.")
    }

    if(is.null(m))
      m <- tmp
    else{
      m <- cbind(m, tmp)
    }
  }

  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]

  norm_method <- match.arg(norm_method)

  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if (requireNamespace("monocle", quietly = TRUE)) {
    if(norm_method == 'vstExprs' && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE) {
      m = monocle::vstExprs(cds, expr_matrix=m)
    }else if(norm_method == 'log') {
      m = log10(m+pseudocount)
    }
  } else {
    warning("Cannot find monocle 'monocle' is not installed.")
  }


  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>scale_max] = scale_max
  m[m<scale_min] = scale_min

  heatmap_matrix <- m

  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1

  if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }else {
    bks <- seq(-3.1,3.1, length.out = length(hmcols))
  }

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ph <- pheatmap::pheatmap(heatmap_matrix,
                             useRaster = T,
                             cluster_cols=FALSE,
                             cluster_rows=T,
                             show_rownames=F,
                             show_colnames=F,
                             clustering_distance_rows=row_dist,
                             clustering_method = hclust_method,
                             cutree_rows=num_clusters,
                             silent=TRUE,
                             filename=NA,
                             breaks=bks,
                             color=hmcols)
  } else {
    warning("Cannot create heatmap. 'pheatmap' is not installed.")
  }


  annotation_col <- data.frame(Branch=factor(rep(rep(branch_label, each = 100))))
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  col_gaps_ind <- c(1:(length(branches) - 1)) * 100

  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(Biobase::fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }


  if (use_gene_short_name == TRUE) {
    if (is.null(Biobase::fData(cds)$gene_short_name) == FALSE) {
      feature_label <- as.character(Biobase::fData(cds)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)

      row_ann_labels <- as.character(Biobase::fData(cds)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }

  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels


  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

  if(!(cluster_rows)) {
    annotation_row <- NA
  }

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ph_res <- pheatmap::pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                                 useRaster = T,
                                 cluster_cols = FALSE,
                                 cluster_rows = cluster_rows,
                                 show_rownames=show_rownames,
                                 show_colnames=F,
                                 #scale="row",
                                 clustering_distance_rows=row_dist, #row_dist
                                 clustering_method = hclust_method, #ward.D2
                                 cutree_rows=num_clusters,
                                 # cutree_cols = 2,
                                 annotation_row=annotation_row,
                                 annotation_col=annotation_col,
                                 gaps_col = col_gaps_ind,
                                 treeheight_row = 20,
                                 breaks=bks,
                                 fontsize = 12,
                                 color=hmcols,
                                 silent=TRUE,
                                 border_color = NA,
                                 filename=NA
    )
  } else {
    warning("Cannot create heatmap. 'pheatmap' is not installed.")
  }



  # ============================================================================
  # prepare data
  wide.res <- cbind(heatmap_matrix,annotation_row) %>%
    data.frame(.,check.names = FALSE) %>%
    dplyr::mutate(gene = rownames(.),.before = 1) %>%
    dplyr::rename(cluster = Cluster)

  # wide to long
  df <- reshape2::melt(wide.res,
                       id.vars = c('cluster','gene'),
                       variable.name = 'cell_type',
                       value.name = 'norm_value') %>%
    dplyr::mutate(cell_type = as.numeric(as.character(cell_type)))

  # add cluster name
  df$cluster_name <- paste('cluster ',df$cluster,sep = '')

  # add gene number
  cl.info <- data.frame(table(wide.res$cluster)) %>%
    dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>%
    dplyr::arrange(Var1)

  id <- unique(df$cluster_name)
  purrr::map_df(seq_along(id),function(x){
    tmp <- df %>%
      dplyr::filter(cluster_name == id[x])

    tmp %>%
      dplyr::mutate(cluster_name = paste(cluster_name," (",cl.info$Freq[x],")",sep = ''))
  }) -> df

  # cluster order
  df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",cl.info$Var1,
                                                           " (",cl.info$Freq,")",sep = ''))

  # return
  prepared_data <- list(wide.res = wide.res,
                        long.res = df,
                        type = "monocle",
                        geneMode = "all",
                        geneType = "branched",
                        pseudotime = annotation_col$Branch)

  if (return_heatmap == TRUE){
    grid::grid.rect(gp=grid::gpar("fill", col=NA))
    grid::grid.draw(ph_res$gtable)
    return(ph_res)
  }else{
    return(prepared_data)
  }
}



#' Calculate and return a smoothed pseudotime matrix for the given gene list
#'
#' This function takes in a monocle3 object and returns a smoothed pseudotime
#' matrix for the
#' given gene list, either in counts or normalized form. The function first matches
#' the gene list
#' with the rownames of the SummarizedExperiment object, and then orders the
#' pseudotime information.
#' The function then uses smooth.spline to apply smoothing to the data. Finally,
#' the function normalizes
#' the data by subtracting the mean and dividing by the standard deviation for each row.
#'
#' @param cds_obj A monocle3 object
#' @param assays Type of assay to be used for the analysis, either "counts" or "normalized"
#' @param gene_list A vector of gene names
#' @return A smoothed pseudotime matrix for the given gene list
#' @importFrom SummarizedExperiment rowData
#' @export
pre_pseudotime_matrix <- function(cds_obj = NULL,
                                  assays = c("counts","normalized"),
                                  gene_list = NULL){
  assays <- match.arg(assays,c("counts","normalized"))

  # # choose assays type
  # if (requireNamespace("monocle3", quietly = TRUE)) {
  #   if(assays == "counts"){
  #     pt.matrix <- monocle3::exprs(cds_obj)[match(gene_list,rownames(SummarizedExperiment::rowData(cds_obj))),
  #                                           order(monocle3::pseudotime(cds_obj))]
  #   }else if(assays == "normalized"){
  #     pt.matrix <- monocle3::normalized_counts(cds_obj, norm_method = "log")[match(gene_list,rownames(SummarizedExperiment::rowData(cds_obj))),
  #                                                                            order(monocle3::pseudotime(cds_obj))]
  #   }
  # } else {
  #   warning("Cannot create monocle3 'monocle3' is not installed.")
  # }

  if(assays == "counts"){
    pt.matrix <- exprs(cds_obj)[match(gene_list,rownames(SummarizedExperiment::rowData(cds_obj))),
                                order(pseudotime(cds_obj))]
  }else if(assays == "normalized"){
    pt.matrix <- normalized_counts(cds_obj, norm_method = "log")[match(gene_list,rownames(SummarizedExperiment::rowData(cds_obj))),
                                                                 order(pseudotime(cds_obj))]
  }

  pt.matrix <- t(apply(pt.matrix,1,function(x){stats::smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- gene_list
  colnames(pt.matrix) <- 1:ncol(pt.matrix)
  return(pt.matrix)
}





#' This is a test data for this package
#' test data describtion
#'
#' @name BEAM_res
#' @docType data
#' @author JunZhang
"BEAM_res"

#' This is a test data for this package
#' test data describtion
#'
#' @name sig_gene_names
#' @docType data
#' @author JunZhang
"sig_gene_names"
