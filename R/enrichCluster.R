#' @name enrichCluster
#' @author JunZhang
#' @title using enrichCluster to do GO/KEGG enrichment analysis for multiple cluster genes
#'
#' @param object clusterData object, default NULL.
#' @param type which type to choose for enrichment "BP","MF","CC" or "KEGG".
#' @param OrgDb the annotation data for enrichment, default NULL.
#' @param id.trans whether perform the ID transformation, default TRUE.
#' @param fromType the input ID type, default "SYMBOL".
#' @param toType the ID type for "bitr" function to transform, default c("ENTREZID").
#' @param readable whether make the enrichmemnt results ID readble, default TRUE.
#' @param organism the organism name for "KEGG" enrichment,mouse("mmu"), human("hsa"), default NULL.
#' @param pvalueCutoff pvalueCutoff for enrichment, default 0.05.
#' @param topn the top enrichment results to extract, length one or same with cluster numbers, default 5.
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
                          id.trans = TRUE,
                          fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          readable = TRUE,
                          # for kegg mouse(mmu)
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
    if(id.trans == TRUE){
      gene.ent <- clusterProfiler::bitr(tmp$gene,
                                        fromType = fromType,
                                        toType = toType,
                                        OrgDb = OrgDb)

      tartget.gene = gene.ent$ENTREZID
    }else{
      tartget.gene = tmp$gene
    }

    # GO enrich
    if(type != "KEGG"){
      set.seed(seed)
      ego <- clusterProfiler::enrichGO(gene          = tartget.gene,
                                       keyType       = toType,
                                       OrgDb         = OrgDb,
                                       ont           = type,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 1,
                                       qvalueCutoff  = 0.2,
                                       readable      = readable)
    }else{
      set.seed(seed)
      ego <- clusterProfiler::enrichKEGG(gene          = tartget.gene,
                                         keyType       = "kegg",
                                         organism      = organism,
                                         universe      = NULL,
                                         pvalueCutoff  = 1,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff  = 0.2)

      # transform gene id
      if(readable == TRUE){
        ego <- clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)
      }else{
        ego <- ego
      }
    }

    # to data.frame
    df <- data.frame(ego) %>%
      dplyr::filter(pvalue < pvalueCutoff) %>%
      dplyr::mutate(group = paste("C",x,sep = '')) %>%
      dplyr::arrange(pvalue)

    # whether save all res
    if(!is.null(topn)){
      # top n selection
      if(length(topn) == 1){
        top = topn
      }else{
        top = topn[x]
      }

      df <- df %>%
        # dplyr::select(group,Description,pvalue) %>%
        dplyr::slice_head(n = top)

      # add gene enrich ratio
      df1 <- purrr::map_df(1:nrow(df),function(x){
        tmp1 <- df[x,]
        size = unlist(strsplit(as.character(tmp1$GeneRatio),split = '/'))
        tmp1$ratio <- (as.numeric(size[1])/as.numeric(size[2]))*100
        return(tmp1)
      })

      # select columns
      df2 <- df1 %>%
        dplyr::select(group,Description,pvalue,ratio)
    }else{
      df2 <- df
    }

    # results
    return(df2)
  })
}

#' This is a test data for this package
#' test data describtion
#'
#' @name net
#' @docType data
#' @author Junjun Lao
"net"
