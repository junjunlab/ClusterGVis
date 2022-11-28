#' @name clusterData
#' @author JunZhang
#' @title using clusterData to cluster genes
#'
#' @param exp the gene normalized expression data(fpkm/rpkm/tpm), default NULL.
#' @param cluster.method which cluster method to choose, mfuzz, kmeans or wgcna, default "mfuzz".
#' details see package Mfuzz. If choose "kmeans", the row_km argument will be used to cluster genes
#' in ComplexHeatmap::Heatmap function.
#' @param object the WGCNA::blockwiseModules returned object, default NULL.
#' @param min.std min std for filter genes, this argument is in mfuzz function, default 0.
#' @param cluster.num the number clusters, default NULL.
#' @param seed set a seed for cluster analysis in mfuzz or Heatmap function, default 5201314.
#' @param scaleData whether do Z-score for expression data, default TRUE.
#'
#' @return clusterData return a list including wide-shape and long-shape clustered results.
#' @export
#'
#' @examples
#'
#' data("exps")
#' # mfuzz
#' cm <- clusterData(exp = exps,
#'                   cluster.method = "mfuzz",
#'                   cluster.num = 8)
#'
#' # kmeans
#' ck <- clusterData(exp = exps,
#'                   cluster.method = "kmeans",
#'                   cluster.num = 8)
#'

globalVariables(c('.', 'cluster', 'cluster2', 'cluster_name','modulecol'))
clusterData <- function(exp = NULL,
                        scaleData = TRUE,
                        cluster.method = c("mfuzz","kmeans","wgcna"),
                        object = NULL,
                        min.std = 0,
                        cluster.num = NULL,
                        seed = 5201314){
  ComplexHeatmap::ht_opt(message = FALSE)

  # choose method
  cluster.method <- match.arg(cluster.method)

  if(cluster.method == "mfuzz"){
    # ========================
    # cluster data
    # myset <- methods::new("ExpressionSet",exprs = as.matrix(exp))
    myset <- Biobase::ExpressionSet(assayData = as.matrix(exp))
    myset <- Mfuzz::filter.std(myset,min.std = min.std,visu = FALSE)

    # whether zsocre data
    if(scaleData == TRUE){
      myset <- Mfuzz::standardise(myset)
    }else{
      myset <- myset
    }

    cluster_number <- cluster.num
    m <- Mfuzz::mestimate(myset)

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
      as.data.frame() %>%
      dplyr::mutate(gene = rownames(.))

    # get gene membership
    lapply(1:cluster.num, function(x){
      ms <- mem %>% dplyr::filter(cluster2 == x)
      res <- data.frame(membership = ms[[x]],gene = ms$gene,cluster2 = ms$cluster2)
    }) %>% do.call('rbind',.) -> membership_info

    # get normalized data
    dnorm <- cbind(myset@assayData$exprs,cluster = mfuzz_res$cluster) %>%
      as.data.frame() %>%
      dplyr::mutate(gene = rownames(.))

    # merge membership info and normalized data
    final_res <- merge(dnorm,membership_info,by = 'gene') %>%
      dplyr::select(-cluster2) %>%
      dplyr::arrange(cluster)

    # wide to long
    df <- reshape2::melt(final_res,
                         id.vars = c('cluster','gene','membership'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')

    # add cluster name
    df$cluster_name <- paste('cluster ',df$cluster,sep = '')

    # cluster order
    df$cluster_name <- factor(df$cluster_name,levels = paste('cluster ',1:cluster.num,sep = ''))

    # add gene number
    cltn <- table(final_res$cluster)
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn],")",sep = ''))
    }) -> df

    # return
    return(list(wide.res = final_res,
                long.res = df,
                type = cluster.method))
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
    set.seed(seed)
    ht = ComplexHeatmap::Heatmap(hclust_matrix,
                                 show_row_names = F,
                                 show_row_dend = F,
                                 row_km = cluster.num)

    # gene order
    ht = ComplexHeatmap::draw(ht)
    row.order = ComplexHeatmap::row_order(ht)

    # get index
    purrr::map_df(1:length(names(row.order)),function(x){
      data.frame(od = row.order[[x]],
                 id = names(row.order)[x])
    }) -> od.res

    cl.info <- data.frame(table(od.res$id))

    # reorder matrix
    m <- hclust_matrix[od.res$od,]

    # add cluster and gene.name
    wide.r <- m %>%
      data.frame() %>%
      dplyr::mutate(gene = rownames(.),
                    cluster = od.res$id) %>%
      dplyr::arrange(cluster)

    # wide to long
    df <- reshape2::melt(wide.r,
                         id.vars = c('cluster','gene'),
                         variable.name = 'cell_type',
                         value.name = 'norm_value')

    # add cluster name
    df$cluster_name <- paste('cluster ',df$cluster,sep = '')

    # cluster order
    df$cluster_name <- factor(df$cluster_name,levels = paste('cluster ',1:cluster.num,sep = ''))

    # add gene number
    cltn <- table(wide.r$cluster)
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn],")",sep = ''))
    }) -> df

    # return
    return(list(wide.res = wide.r,
                long.res = df,
                type = cluster.method))
  }else if(cluster.method == "wgcna"){
    # =====================================
    net <- object
    cinfo <- data.frame(cluster = net$colors + 1,
                        modulecol = WGCNA::labels2colors(net$colors))

    expm <- data.frame(t(scale(exp)))
    expm$gene <- rownames(expm)
    final.res <- data.frame(cbind(expm,cinfo)) %>%
      dplyr::arrange(cluster)

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
    df$cluster_name <- factor(df$cluster_name,levels = paste('cluster ',1:nrow(cl.info),sep = ''))

    # add gene number
    cltn <- table(final.res$cluster)
    purrr::map_df(unique(df$cluster_name),function(x){
      tmp <- df %>%
        dplyr::filter(cluster_name == x)

      cn = as.numeric(unlist(strsplit(as.character(x),split = "cluster "))[2])

      tmp %>%
        dplyr::mutate(cluster_name = paste(cluster_name," (",cltn[cn]," ",unique(tmp$modulecol),")",sep = ''))
    }) -> df

    # return
    return(list(wide.res = final.res,
                long.res = df,
                type = cluster.method))
  }else{
    print("supply with mfuzz, kmeans or wgcna !")
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
