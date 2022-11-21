#' @name enrichCluster
#' @author JunZhang
#' @title using enrichCluster to do GO/KEGG enrichment analysis for multiple cluster genes
#'
#' @param object clusterData object, default NULL.
#' @param type which type to choose for enrichment "BP","MF","CC" or "KEGG".
#' @param OrgDb the annotation data for enrichment, default NULL.
#' @param organism the organism name for kegg enrichment, default NULL.
#' @param pvalueCutoff pvalueCutoff for enrichment, default 0.05.
#' @param topn the top enrichment results to extract, default 5.
#' @param seed the enrichment seed, default 5201314.
#'
#' @import org.Mm.eg.db
#'
#' @return a data.frame.
#' @export
#'
#' @examples
#' \dontrun{
#' library(org.Mm.eg.db)
#'
#' data("exps")
#'
#' # mfuzz
#' ck <- clusterData(exp = exps,
#'                   cluster.method = "kmeans",
#'                   cluster.num = 8)
#'
#' # enrich for clusters
#' enrich <- enrichCluster(object = ck,
#'                         OrgDb = org.Mm.eg.db,
#'                         type = "BP",
#'                         topn = 5)
#' }
globalVariables(c('Description', 'group', 'pvalue'))
enrichCluster <- function(object = NULL,
                          type = c("BP","MF","CC","KEGG"),
                          OrgDb = NULL,
                          # for kegg
                          organism = "hsa",
                          pvalueCutoff  = 0.05,
                          topn = 5,
                          seed = 5201314){
  type <- match.arg(type)

  # get data
  enrich.data <- object$wide.res

  # loop for enrich
  purrr::map_df(1:length(unique(enrich.data$cluster)),function(x){
    # filter
    tmp <- enrich.data %>%
      dplyr::filter(cluster == x)

    # =============================================
    # enrich

    # entriz id transformation
    gene.ent <- clusterProfiler::bitr(tmp$gene,
                                      fromType = "SYMBOL",
                                      toType = c("ENTREZID"),
                                      OrgDb = OrgDb)

    # GO enrich
    if(type != "KEGG"){
      set.seed(seed)
      ego <- clusterProfiler::enrichGO(gene          = gene.ent$ENTREZID,
                                       keyType       = "ENTREZID",
                                       OrgDb         = OrgDb,
                                       ont           = type,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = pvalueCutoff,
                                       qvalueCutoff  = 0.2,
                                       readable      = TRUE)
    }else{
      set.seed(seed)
      ego <- clusterProfiler::enrichKEGG(gene          = gene.ent$ENTREZID,
                                         keyType       = "kegg",
                                         organism      = organism,
                                         pvalueCutoff  = pvalueCutoff,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff  = 0.2)

      # transform gene id
      ego <- clusterProfiler::setReadable(ego,OrgDb = OrgDb)
    }

    # to data.frame
    df <- data.frame(ego) %>%
      dplyr::mutate(group = paste("C",x,sep = '')) %>%
      dplyr::arrange(pvalue)

    # whether save all res
    if(!is.null(topn)){
      df <- df %>%
        dplyr::select(group,Description,pvalue) %>%
        dplyr::slice_head(n = topn)
    }else{
      df <- df
    }

    # results
    return(df)
  })
}
