globalVariables(c('.', 'cluster', 'cluster2', 'cluster_name','modulecol',"cl.x", "cl.y"))

#' Cluster Data Based on Different Methods
#' @author JunZhang
#'
#' This function performs clustering on input data using one of four methods:
#' **mfuzz**, **TCseq**, **kmeans**, or **wgcna**. The clustering results include
#' metadata, normalized data, and cluster memberships.
#'
#' @param obj An input object that can take one of two types:
#'        - A **cell_data_set** object for trajectory analysis.
#'        - A **matrix** or **data.frame** containing expression data.
#' @param scaleData Logical. Whether to scale the data (e.g., z-score normalization).
#' @param cluster.method Character. Clustering method to use.
#'        Options are one of `"mfuzz"`, `"TCseq"`, `"kmeans"`, or `"wgcna"`.
#' @param TCseq_params_list A list of additional parameters passed to the `TCseq::timeclust` function.
#' @param object A pre-calculated object required when using `"wgcna"` as the clustering method.
#' @param min.std Numeric. Minimum standard deviation for filtering expression data.
#' @param cluster.num Integer. The number of clusters to identify.
#' @param subcluster A numeric vector of specific cluster IDs to include in the results.
#'        If `NULL`, all clusters are included.
#' @param seed An integer seed for reproducibility in clustering operations.
#' @param ... Additional arguments passed to internal functions such as `pre_pseudotime_matrix`.
#'
#' @return A list containing the following clustering results:
#'   - **wide.res**: A wide-format data frame with clusters and normalized expression levels.
#'   - **long.res**: A long-format data frame for visualizations, containing cluster information,
#'                   normalized values, cluster names, and memberships.
#'   - **cluster.list**: A list where each element contains genes belonging to a specific cluster.
#'   - **type**: The clustering method used (`"mfuzz"`, `"TCseq"`, `"kmeans"`, or `"wgcna"`).
#'   - **geneMode**: Currently set to `"none"` (reserved for future use).
#'   - **geneType**: Currently set to `"none"` (reserved for future use).
#'
#' @details
#' Depending on the selected `cluster.method`, different clustering algorithms are used:
#'   - **`"mfuzz"`**: Applies Mfuzz soft clustering method, suitable for identifying overlapping clusters.
#'   - **`"TCseq"`**: Uses TCseq clustering for time-series expression data with support for additional parameters.
#'   - **`"kmeans"`**: Employs standard k-means clustering via base R's `stats::kmeans`.
#'   - **`"wgcna"`**: Leverages pre-calculated WGCNA (Weighted Gene Co-expression Network Analysis) networks.
#'
#' The function is designed to be flexible, allowing preprocessing (e.g., filtering by `min.std`),
#' scaling the data (`scaleData = TRUE`), and generating results compatible with data visualization pipelines.
#'
#' @section WGCNA Clustering:
#' If the **WGCNA** method is selected, the `object` parameter must contain a pre-calculated WGCNA network object.
#' This is typically obtained using the WGCNA package functions.
#'
#' @section Subsetting Clusters:
#' Use the `subcluster` parameter to focus on specific clusters. Cluster IDs not included in the
#' `subcluster` vector will be excluded from the final results.
#'
#' @importFrom utils modifyList
#' @importFrom stats kmeans
#' @importFrom Biobase exprs
#'
#' @export
#'
#' @examples
#'
#' data("exps")
#'
#' # kmeans
#' ck <- clusterData(obj = exps,
#'                   cluster.method = "kmeans",
#'                   cluster.num = 8)
#'
clusterData <- function(obj = NULL,
                        scaleData = TRUE,
                        cluster.method = c("mfuzz","TCseq","kmeans","wgcna"),
                        TCseq_params_list = list(),
                        object = NULL,
                        min.std = 0,
                        cluster.num = NULL,
                        subcluster = NULL,
                        seed = 5201314,...){
  ComplexHeatmap::ht_opt(message = FALSE)

  # check datatype
  cls <- class(obj)
  # pkg <- attr(cls,"package")

  if("cell_data_set" %in% cls){
    extra_params <- list(cds_obj = obj,assays = "counts",...)

    exp <- do.call(pre_pseudotime_matrix,extra_params)
  }else if("matrix" %in% cls | "data.frame" %in% cls){
    exp <- obj
  }

  # choose method
  cluster.method <- match.arg(cluster.method)

  # check clusting method
  if(cluster.method %in% c("mfuzz","TCseq")){

    if(cluster.method == "mfuzz"){
      # ==========================================================================
      # mfuzz
      # ==========================================================================
      # cluster data
      # myset <- methods::new("ExpressionSet",exprs = as.matrix(exp))
      myset <- Biobase::ExpressionSet(assayData = as.matrix(exp))
      myset <- filter.std(myset,min.std = min.std,visu = FALSE)

      # whether zsocre data
      if(scaleData == TRUE){
        myset <- standardise(myset)
      }else{
        myset <- myset
      }

      cluster_number <- cluster.num
      m <- mestimate(myset)

      # cluster step
      set.seed(seed)
      # mfuzz_res <- Mfuzz::mfuzz(myset, c = cluster_number, m = m)

      # user define mfuzz
      mfuzz1 <- function(eset,centers,m,...){
        cl <- e1071::cmeans(Biobase::exprs(eset),centers = centers,method = "cmeans",m = m,...)
      }

      mfuzz_res <- mfuzz1(myset, c = cluster_number, m = m)
      # ========================
      # get clustered data
      mtx <- Biobase::assayData(myset)
      mtx <- mtx$exprs
      raw_cluster_anno <- cbind(mtx,cluster = mfuzz_res$cluster)

      # membership
      mem <- cbind(mfuzz_res$membership,cluster2 = mfuzz_res$cluster) %>%
        data.frame(check.names = FALSE) %>%
        dplyr::mutate(gene = rownames(.))

      # get gene membership
      lapply(1:cluster.num, function(x){
        ms <- mem %>% dplyr::filter(cluster2 == x)
        res <- data.frame(membership = ms[[x]],gene = ms$gene,cluster2 = ms$cluster2,
                          check.names = FALSE)
      }) %>% do.call('rbind',.) -> membership_info

      # get normalized data
      dnorm <- cbind(myset@assayData$exprs,cluster = mfuzz_res$cluster) %>%
        data.frame(check.names = FALSE) %>%
        dplyr::mutate(gene = rownames(.))

      # merge membership info and normalized data
      final_res <- merge(dnorm,membership_info,by = 'gene') %>%
        dplyr::select(-cluster2) %>%
        dplyr::arrange(cluster)
    }else if(cluster.method == "TCseq"){
      # ==========================================================================
      # TCseq
      # ==========================================================================
      exp <- filter.std(exp,min.std = min.std,visu = FALSE)

      # tca <- TCseq::timeclust(x = as.matrix(exp), algo = "cm",
      #                         k = cluster.num, standardize = scaleData)

      tca <- do.call(TCseq::timeclust,
                     modifyList(list(x = as.matrix(exp),
                                     algo = "cm",
                                     k = cluster.num,
                                     standardize = scaleData),
                                TCseq_params_list))

      dt <- data.frame(tca@data) %>%
        tibble::rownames_to_column(var = "gene")

      cluster <- data.frame(cl = tca@cluster,check.names = F) %>%
        tibble::rownames_to_column(var = "gene")

      membership <- data.frame(tca@membership,check.names = F) %>%
        tibble::rownames_to_column(var = "gene") %>%
        reshape2::melt(id.vars = "gene",variable.name = "cl",value.name = "membership")

      anno <- cluster %>%
        dplyr::left_join(y = membership,by = "gene") %>%
        dplyr::filter(cl.x == cl.y) %>%
        dplyr::select(-cl.y) %>%
        dplyr::rename(cluster = cl.x)

      final_res <- dt %>%
        dplyr::left_join(y = anno,by = "gene")

    }


    # ==========================================================================
    # whether subset clusters
    if(!is.null(subcluster)){
      final_res <- final_res %>% dplyr::filter(cluster %in% subcluster)
    }

    # wide to long
    df <- reshape2::melt(final_res,
                         id.vars = c('cluster','gene','membership'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')

    # add cluster name
    df$cluster_name <- paste('cluster ',df$cluster,sep = '')

    # add gene number
    cltn <- table(final_res$cluster)
    cl.info <- data.frame(table(final_res$cluster))
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn],")",sep = ''))
    }) -> df

    # cluster order
    df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info),
                                                             " (",cl.info$Freq,")",sep = ''))

    # cluster list
    wide <- final_res
    wide$cluster <- paste("C",wide$cluster,sep = "")
    cluster.list <- split(wide$gene,wide$cluster)

    # return
    return(list(wide.res = final_res,
                long.res = df,
                cluster.list = cluster.list,
                type = cluster.method,
                geneMode = "none",
                geneType = "none"))
  }else if(cluster.method == "kmeans"){
    # ==========================================================================
    # using complexheatmap cluster genes
    exp <- filter.std(exp,min.std = min.std,visu = FALSE)

    # whether zsocre data
    if(scaleData == TRUE){
      hclust_matrix <- exp %>% t() %>% scale() %>% t()
    }else{
      hclust_matrix <- exp
    }

    # plot
    # set.seed(seed)
    # ht <- ComplexHeatmap::Heatmap(hclust_matrix,
    #                               show_row_names = F,
    #                               show_row_dend = F,
    #                               show_column_names = F,
    #                               row_km = cluster.num)
    #
    # # gene order
    # ht = ComplexHeatmap::draw(ht)
    # row.order = ComplexHeatmap::row_order(ht)
    #
    # # get index
    # if(is.null(cluster.num) | cluster.num == 1){
    #   od.res <- data.frame(od = row.order,id = 1)
    # }else{
    #   # get index
    #   purrr::map_df(1:length(names(row.order)),function(x){
    #     data.frame(od = row.order[[x]],
    #                id = as.numeric(names(row.order)[x]),
    #                check.names = FALSE)
    #   }) -> od.res
    # }

    # add kmeans func n stats
    km <- stats::kmeans(x = hclust_matrix,centers = cluster.num,nstart = 10)

    od.res <- data.frame(od = match(names(km$cluster),rownames(hclust_matrix)),
                         id = as.numeric(km$cluster),
                         check.names = FALSE)


    cl.info <- data.frame(table(od.res$id),check.names = FALSE)

    # reorder matrix
    m <- hclust_matrix[od.res$od,]

    # add cluster and gene.name
    wide.r <- m %>%
      data.frame(check.names = FALSE) %>%
      dplyr::mutate(gene = rownames(.),
                    cluster = od.res$id) %>%
      dplyr::arrange(cluster)

    # whether subset clusters
    if(!is.null(subcluster)){
      wide.r <- wide.r %>% dplyr::filter(cluster %in% subcluster)
    }

    # wide to long
    df <- reshape2::melt(wide.r,
                         id.vars = c('cluster','gene'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')

    # add cluster name
    df$cluster_name <- paste('cluster ',df$cluster,sep = '')

    # add gene number
    cltn <- table(wide.r$cluster)
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn],")",sep = ''))
    }) -> df

    # cluster order
    df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info),
                                                             " (",cl.info$Freq,")",sep = ''))
    # cluster list
    wide <- wide.r
    wide$cluster <- paste("C",wide$cluster,sep = "")
    cluster.list <- split(wide$gene,wide$cluster)

    # return
    return(list(wide.res = wide.r,
                long.res = df,
                cluster.list = cluster.list,
                type = cluster.method,
                geneMode = "none",
                geneType = "none"))
  }else if(cluster.method == "wgcna"){
    # =====================================
    net <- object
    cinfo <- data.frame(cluster = net$colors + 1,
                        modulecol = WGCNA::labels2colors(net$colors),
                        check.names = FALSE)

    expm <- data.frame(t(scale(exp)))
    expm$gene <- rownames(expm)
    final.res <- data.frame(cbind(expm,cinfo)) %>%
      dplyr::arrange(cluster)

    # whether subset clusters
    if(!is.null(subcluster)){
      final.res <- final.res %>% dplyr::filter(cluster %in% subcluster)
    }

    # =====================================
    cl.info <- data.frame(table(final.res$cluster))

    # wide to long
    df <- reshape2::melt(final.res,
                         id.vars = c('cluster','gene','modulecol'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')

    # add cluster name
    df$cluster_name <- paste('cluster ',df$cluster,sep = '')

    # cluster order
    # df$cluster_name <- factor(df$cluster_name,levels = paste('cluster ',1:nrow(cl.info),sep = ''))

    # add gene number
    cltn <- table(final.res$cluster)
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn]," ",unique(tmp$modulecol),")",sep = ''))
    }) -> df

    # cluster order
    # df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info),
    #                                                          " (",cl.info$Freq,")",sep = ''))

    df$cluster_name <- factor(df$cluster_name,levels = paste("cluster ",1:nrow(cl.info), " (",cl.info$Freq," ",unique(df$modulecol),")",sep = ''))

    # cluster list
    wide <- final.res
    wide$cluster <- paste("C",wide$cluster,sep = "")
    cluster.list <- split(wide$gene,wide$cluster)

    # return
    return(list(wide.res = final.res,
                long.res = df,
                cluster.list = cluster.list,
                type = cluster.method,
                geneMode = "none",
                geneType = "none"))
  }else{
    message("supply with mfuzz, kmeans or wgcna !")
  }
}


###############################
#' This is a test data for this package
#' test data describtion
#'
#' @name exps
#' @docType data
#' @author Junjun Lao
"exps"
