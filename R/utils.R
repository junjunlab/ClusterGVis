globalVariables(c("newCellDataSet"))


# from mfuzz
standardise <- function(eset) {
  data <- Biobase::exprs(eset)
  for (i in seq_len(dim(data)[[1]])) {
    data[i, ] <-
      (data[i, ] - base::mean(data[i, ], na.rm = TRUE)) / stats::sd(data[i, ],
                                                                    na.rm =
                                                                      TRUE
      )
  }
  Biobase::exprs(eset) <- data
  eset
}


# from mfuzz
mestimate <- function(eset) {
  N <- dim(Biobase::exprs(eset))[[1]]
  D <- dim(Biobase::exprs(eset))[[2]]
  m.sj <-
    1 + (1418 / N + 22.05) * D^(-2) + (12.33 / N + 0.243) * D^(-0.0406 *
                                                                 log(N) - 0.1134)
  return(m.sj)
}


# Test whether a matrix is one of our supported sparse matrices
# author https://github.com/cole-trapnell-lab/monocle3
is_sparse_matrix <- function(x) {
  class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix")
}



size_factors <- function(cds) {
  stopifnot(methods::is(cds, "cell_data_set"))

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required. Please install it.")
  }

  sf <- SingleCellExperiment::colData(cds)$Size_Factor
  names(sf) <- colnames(SingleCellExperiment::counts(cds))
  sf
}




# ==============================================================================

#' Cell Data Set Class
#' @name cell_data_set-class
#' @rdname cell_data_set-class
#' @exportClass cell_data_set
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setClass(
  "cell_data_set",
  contains = "SingleCellExperiment",
  slots = c(
    reduce_dim_aux = "SimpleList",
    principal_graph_aux = "SimpleList",
    principal_graph = "SimpleList",
    clusters = "SimpleList"
  )
)



# ==============================================================================

#' Generic to extract pseudotime from CDS object
#'
#' @author https://github.com/cole-trapnell-lab/monocle3
#'
#'
#' @param x A cell_data_set object.
#' @param reduction_method Reduced dimension to extract pseudotime for.
#'
#'
#' @return A named numeric vector containing pseudotime values for each cell.
#'   The vector has the following properties:
#'   \itemize{
#'     \item Names correspond to cell identifiers (column names from the
#'     expression matrix)
#'     \item Values are numeric pseudotime values representing trajectory
#'     progress
#'     \item Values typically range from 0 (trajectory start) to maximum
#'     trajectory length
#'     \item Higher values indicate cells further along the developmental
#'      trajectory
#'   }
#'
#'
#'
setGeneric("pseudotime", function(x, reduction_method = "UMAP") {
  standardGeneric("pseudotime")
})



#' @rdname pseudotime
#'
setMethod("pseudotime", "cell_data_set",
          function(x, reduction_method = "UMAP") {
            value <- x@principal_graph_aux[[reduction_method]]$pseudotime[
              colnames(exprs(x))]
            if (is.null(value)) {
              stop(
                "No pseudotime calculated for reduction_method = ",
                reduction_method,
                ". Please first run ",
                "order_cells with reduction_method = ",
                reduction_method,
                "."
              )
            }
            return(value)
          })



#' Generic to access cds count matrix
#'
#' @author https://github.com/cole-trapnell-lab/monocle3
#'
#' @param x A cell_data_set object.
#'
#'
#' @return Count matrix.
#'
#'
setGeneric("exprs", function(x) {
  standardGeneric("exprs")
})



#' Method to access cds count matrix
#' @param x A cell_data_set object.
#'
#' @return Count matrix.
#'
#'
setMethod("exprs", "cell_data_set", function(x) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required. Please install it.")
  }

  value <- SummarizedExperiment::assays(x)$counts
  return(value)
})


# ==============================================================================

# Return a size-factor normalized and (optionally) log-transformed expression
#
# @author https://github.com/cole-trapnell-lab/monocle3
#
normalized_counts <- function(cds,
                              norm_method = c("log", "binary", "size_only"),
                              pseudocount = 1) {
  norm_method <- match.arg(norm_method,
                           choices = c("log", "binary", "size_only"))

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required. Please install it.")
  }

  norm_mat <- SingleCellExperiment::counts(cds)
  if (norm_method == "binary") {
    # The '+ 0' coerces the matrix to type numeric. It's possible
    # to use 'as.numeric(norm_mat > 0)' but the matrix
    # attributes disappear...
    norm_mat <- (norm_mat > 0) + 0
    if (is_sparse_matrix(norm_mat)) {
      norm_mat <- methods::as(norm_mat, "dgCMatrix")
    }
  } else {
    if (is_sparse_matrix(norm_mat)) {
      norm_mat@x <- norm_mat@x / rep.int(size_factors(cds), diff(norm_mat@p))
      if (norm_method == "log") {
        if (pseudocount == 1) {
          norm_mat@x <- log10(norm_mat@x + pseudocount)
        } else {
          stop("Pseudocount must equal 1 with sparse expression matrices")
        }
      }
    } else {
      norm_mat <- Matrix::t(Matrix::t(norm_mat) / size_factors(cds))
      if (norm_method == "log") {
        norm_mat@x <- log10(norm_mat + pseudocount)
      }
    }
  }
  return(norm_mat)
}


minSpanningTree <- function( cds ) {
  stopifnot(methods::is( cds, "CellDataSet" ) )
  cds@minSpanningTree
}

`minSpanningTree<-` <- function( cds, value ) {
  stopifnot(methods::is( cds, "CellDataSet" ) )
  cds@minSpanningTree <- value
  methods::validObject( cds )
  cds
}

buildBranchCellDataSet2 <- function(cds,
                                    progenitor_method = c('sequential_split', 'duplicate'),
                                    branch_states = NULL,
                                    branch_point = 1,
                                    branch_labels = NULL,
                                    stretch = TRUE){
  # TODO: check that branches are on the same paths
  if(is.null(Biobase::pData(cds)$State) | is.null(Biobase::pData(cds)$Pseudotime))
    stop('Please first order the cells in pseudotime using orderCells()')
  if(is.null(branch_point) & is.null(branch_states))
    stop('Please either specify the branch_point or branch_states to select subset of cells')
  #if(ncol(cds@reducedDimS) != ncol(cds))
  #  stop('You probably used clusterCells function which should be used together with buildBranchCellDataSet, try re-run reduceDimension without clustering cells again')

  if (!is.null(branch_labels) & !is.null(branch_states)) {
    if(length(branch_labels) != length(branch_states))
      stop("length of branch_labels doesn't match with that of branch_states")
    branch_map <- stats::setNames(branch_labels, as.character(branch_states))
  }

  if(cds@dim_reduce_type == "DDRTree") {
    pr_graph_cell_proj_mst <- minSpanningTree(cds)
  }
  else {
    pr_graph_cell_proj_mst <- cds@auxOrderingData[[cds@dim_reduce_type]]$cell_ordering_tree
  }

  root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root_state <- Biobase::pData(cds)[root_cell,]$State
  #root_state <- V(pr_graph_cell_proj_mst)[root_cell,]$State

  pr_graph_root <- subset(Biobase::pData(cds), State == root_state)

  if (cds@dim_reduce_type == "DDRTree"){
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
    root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root),]
  }else{
    root_cell_point_in_Y <- row.names(pr_graph_root)
  }

  root_cell <- names(which(igraph::degree(pr_graph_cell_proj_mst,
                                          v = root_cell_point_in_Y,
                                          mode = "all")==1, useNames = TRUE))[1]

  paths_to_root <- list()
  if (is.null(branch_states) == FALSE){

    # If the user didn't specify a branch point,
    # let's walk back from the branch states
    for (leaf_state in branch_states){

      curr_cell <- subset(Biobase::pData(cds), State == leaf_state)
      #Get all the nearest cells in Y for curr_cells:

      if (cds@dim_reduce_type == "DDRTree"){
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        curr_cell_point_in_Y <- closest_vertex[row.names(curr_cell),]
      }else{
        curr_cell_point_in_Y <- row.names(curr_cell)
      }

      # Narrow down to a single tip cell in Y:
      curr_cell <- names(which(igraph::degree(pr_graph_cell_proj_mst,
                                              v = curr_cell_point_in_Y,
                                              mode = "all")==1, useNames = TRUE))[1]

      path_to_ancestor <- igraph::shortest_paths(pr_graph_cell_proj_mst,
                                                 curr_cell, root_cell)
      path_to_ancestor <- names(unlist(path_to_ancestor$vpath))

      if (cds@dim_reduce_type == "DDRTree"){
        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        ancestor_cells_for_branch <- row.names(closest_vertex)[which(igraph::V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_ancestor)]
      }else if (cds@dim_reduce_type == "ICA"){
        ancestor_cells_for_branch <- path_to_ancestor
      }
      ancestor_cells_for_branch <- intersect(ancestor_cells_for_branch, colnames(cds))
      paths_to_root[[as.character(leaf_state)]] <- ancestor_cells_for_branch
    }
  }else{
    if(cds@dim_reduce_type == "DDRTree")
      pr_graph_cell_proj_mst <- minSpanningTree(cds)
    else
      pr_graph_cell_proj_mst <- cds@auxOrderingData$ICA$cell_ordering_tree

    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_cell <- mst_branch_nodes[branch_point]
    mst_no_branch_point <- pr_graph_cell_proj_mst - igraph::V(pr_graph_cell_proj_mst)[branch_cell]

    path_to_ancestor <- igraph::shortest_paths(pr_graph_cell_proj_mst, branch_cell, root_cell)
    path_to_ancestor <- names(unlist(path_to_ancestor$vpath))

    #post_branch_cells <- c()
    for (backbone_nei in igraph::V(pr_graph_cell_proj_mst)[suppressWarnings(.nei(branch_cell))]$name){
      descendents <- igraph::bfs(mst_no_branch_point, igraph::V(mst_no_branch_point)[backbone_nei], unreachable=FALSE)
      descendents <- descendents$order[!is.na(descendents$order)]
      descendents <- igraph::V(mst_no_branch_point)[descendents]$name
      if (root_cell %in% descendents == FALSE){
        path_to_root <- unique(c(path_to_ancestor, branch_cell, descendents))

        if (cds@dim_reduce_type == "DDRTree"){
          closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
          path_to_root <- row.names(closest_vertex)[
            which(igraph::V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_root)]
        }else{
          path_to_root <- path_to_root
        }

        closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
        #branch_state <- unique(pData(cds)[backbone_nei, ]$State)[1]

        path_to_root <- intersect(path_to_root, colnames(cds))
        paths_to_root[[backbone_nei]] <- path_to_root
        #post_branch_cells <- c(post_branch_cells, backbone_nei)
      }
    }
  }
  all_cells_in_subset <- c()

  if (is.null(branch_labels) == FALSE){
    if (length(branch_labels) != 2)
      stop("Error: branch_labels must have exactly two entries")
    names(paths_to_root) <- branch_labels
  }

  for (path_to_ancestor in paths_to_root){
    if (length(path_to_ancestor) == 0){
      stop("Error: common ancestors between selected State values on path to root State")
    }
    all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
  }
  all_cells_in_subset <- unique(all_cells_in_subset)

  common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
  # if (length(paths_to_root) > 2){
  #   for (i in seq(3,length(paths_to_root))){
  #     common_ancestor_cells <- intersect(intersect(paths_to_root[i], paths_to_root[i-1]), common_ancestor_cells)
  #   }
  # }

  #when n-center used, this creates problems
  cds <- cds[, row.names(Biobase::pData(cds[,all_cells_in_subset]))] #or just union(ancestor_cells, branch_cells)

  #State <- pData(cds)$State
  Pseudotime <- Biobase::pData(cds)$Pseudotime

  pData <- Biobase::pData(cds)

  if(stretch) {
    max_pseudotime <- -1
    for (path_to_ancestor in paths_to_root){
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime)
      if (max_pseudotime < max_pseudotime_on_path){
        max_pseudotime <- max_pseudotime_on_path
      }
    }

    branch_pseudotime <- max(pData[common_ancestor_cells,]$Pseudotime)
    #ancestor_scaling_factor <- branch_pseudotime / max_pseudotime

    for (path_to_ancestor in paths_to_root){
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime)
      path_scaling_factor <-(max_pseudotime - branch_pseudotime) / (max_pseudotime_on_path - branch_pseudotime)
      if (is.finite(path_scaling_factor)){
        branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
        pData[branch_cells,]$Pseudotime <- ((pData[branch_cells,]$Pseudotime - branch_pseudotime) * path_scaling_factor + branch_pseudotime)
      }
    }
    #pData[common_ancestor_cells,]$Pseudotime <- pData[common_ancestor_cells,]$Pseudotime / max_pseudotime

    pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime
  }
  pData$original_cell_id <- row.names(pData)

  pData$original_cell_id <- row.names(pData)

  if(length(paths_to_root) != 2)
    stop('more than 2 branch states are used!')

  pData[common_ancestor_cells, "Branch"] <- names(paths_to_root)[1] #set progenitors to the branch 1

  progenitor_pseudotime_order <- order(pData[common_ancestor_cells, 'Pseudotime'])

  # if (progenitor_method == 'duplicate') {
  if('duplicate' %in% progenitor_method){
    ancestor_exprs <- Biobase::exprs(cds)[,common_ancestor_cells]
    expr_blocks <- list()

    # Duplicate the expression data
    for (i in seq_len(length(paths_to_root))) { #duplicate progenitors for multiple branches
      if (nrow(ancestor_exprs) == 1)
        exprs_data <- t(as.matrix(ancestor_exprs))
      else exprs_data <- ancestor_exprs

      colnames(exprs_data) <- paste('duplicate', i, seq_len(length(common_ancestor_cells)), sep = '_')
      expr_lineage_data <- Biobase::exprs(cds)[,setdiff(paths_to_root[[i]], common_ancestor_cells)]
      exprs_data <- cbind(exprs_data, expr_lineage_data)
      expr_blocks[[i]] <- exprs_data
    }

    # Make a bunch of copies of the pData entries from the common ancestors
    ancestor_pData_block <- pData[common_ancestor_cells,]

    pData_blocks <- list()

    weight_vec <- c()
    for (i in seq_len(length(paths_to_root))) {
      weight_vec <- c(weight_vec, rep(1, length(common_ancestor_cells)))
      weight_vec_block <- rep(1, length(common_ancestor_cells))

      #pData <- rbind(pData, pData[common_ancestor_cells, ])
      new_pData_block <- ancestor_pData_block
      # new_pData_block$Lineage <- lineage_states[i]
      # new_pData_block$State <- lineage_states[i]

      row.names(new_pData_block) <- paste('duplicate', i, seq_len(length(common_ancestor_cells)), sep = '_')

      pData_lineage_cells <- pData[setdiff(paths_to_root[[i]], common_ancestor_cells),]
      # pData_lineage_cells$Lineage <- lineage_states[i]
      # pData_lineage_cells$State <- lineage_states[i]

      weight_vec_block <- c(weight_vec_block, rep(1, nrow(pData_lineage_cells)))

      weight_vec <- c(weight_vec, weight_vec_block)

      new_pData_block <- rbind(new_pData_block, pData_lineage_cells)
      new_pData_block$Branch <- names(paths_to_root)[i]
      pData_blocks[[i]] <- new_pData_block
    }
    pData <- do.call(rbind, pData_blocks)
    exprs_data <- do.call(cbind, expr_blocks)
  }
  else if(progenitor_method == 'sequential_split') {
    pData$Branch <- names(paths_to_root)[1]

    branchA <- progenitor_pseudotime_order[seq(1, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[branchA], 'Branch'] <- names(paths_to_root)[1]
    branchB <- progenitor_pseudotime_order[seq(2, length(common_ancestor_cells), by = 2)]
    pData[common_ancestor_cells[branchB], 'Branch'] <- names(paths_to_root)[2]

    # Duplicate the root cell to make sure both regression start at pseudotime zero:
    zero_pseudotime_root_cell <- common_ancestor_cells[progenitor_pseudotime_order[1]]
    exprs_data <- cbind(Biobase::exprs(cds), 'duplicate_root' = Biobase::exprs(cds)[, zero_pseudotime_root_cell])
    pData <- rbind(pData, pData[zero_pseudotime_root_cell, ])
    row.names(pData)[nrow(pData)] <- 'duplicate_root'
    pData[nrow(pData), 'Branch'] <- names(paths_to_root)[2]

    weight_vec <- rep(1, nrow(pData))

    for (i in seq_len(length(paths_to_root))){
      path_to_ancestor <- paths_to_root[[i]]
      branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
      pData[branch_cells,]$Branch <- names(paths_to_root)[i]
    }
  }

  pData$Branch <- as.factor(pData$Branch)

  pData$State <- factor(pData$State)
  Size_Factor <- pData$Size_Factor

  fData <- Biobase::fData(cds)

  colnames(exprs_data) <- row.names(pData) #check this
  cds_subset <- monocle::newCellDataSet(as.matrix(exprs_data),
                                        phenoData = methods::new("AnnotatedDataFrame", data = pData),
                                        featureData = methods::new("AnnotatedDataFrame", data = fData),
                                        expressionFamily= cds@expressionFamily,
                                        lowerDetectionLimit = cds@lowerDetectionLimit)
  Biobase::pData(cds_subset)$State <- as.factor(Biobase::pData(cds_subset)$State)
  Biobase::pData(cds_subset)$Size_Factor <- Size_Factor

  cds_subset@dispFitInfo <- cds@dispFitInfo

  return (cds_subset)
}
