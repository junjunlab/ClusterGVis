globalVariables(c('Description', 'group', 'pvalue',"geneID"))

#' @name enrichCluster
#' @author JunZhang
#' @title using enrichCluster to do GO/KEGG enrichment analysis for multiple cluster genes
#'
#' @param object clusterData object, default NULL.
#' @param type which type to choose for enrichment "BP","MF","CC","KEGG" or "ownSet".
#' @param TERM2GENE user input annotation of TERM TO GENE mapping, a data.frame
#' of 2 column with term and gene, default NULL.
#' @param TERM2NAME user input of TERM TO NAME mapping, a data.frame of 2 column
#' with term and name, default NULL.
#' @param OrgDb the annotation data for enrichment, default NULL.
#' @param id.trans whether perform the ID transformation, default TRUE.
#' @param fromType the input ID type, default "SYMBOL".
#' @param toType the ID type for "bitr" function to transform, default c("ENTREZID").
#' @param readable whether make the enrichmemnt results ID readble, default TRUE.
#' @param organism the organism name for "KEGG" enrichment,mouse("mmu"), human("hsa"), default NULL.
#' @param pvalueCutoff pvalueCutoff for enrichment, default 0.05.
#' @param topn the top enrichment results to extract, length one or same with cluster numbers, default 5.
#' @param seed the enrichment seed, default 5201314.
#' @param add.gene whether return genes for output, default FALSE.
#'
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
enrichCluster <- function(object = NULL,
                          type = c("BP","MF","CC","KEGG","ownSet"),
                          TERM2GENE = NULL,
                          TERM2NAME = NULL,
                          OrgDb = NULL,
                          id.trans = TRUE,
                          fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          readable = TRUE,
                          # for kegg mouse(mmu)
                          organism = "hsa",
                          pvalueCutoff  = 0.05,
                          topn = 5,
                          seed = 5201314,
                          add.gene = FALSE){
  type <- match.arg(type)

  # get data
  enrich.data <- object$wide.res

  # check geneType
  if(startsWith(object$geneType,"nounique")){
    split = unlist(strsplit(object$geneType,split = "\\|"))[2]
    enrich.data$gene <- sapply(strsplit(as.character(enrich.data$gene),split = split),"[",1)
  }

  # loop for enrich
  purrr::map_df(1:length(unique(enrich.data$cluster)),function(x){
    # filter
    tmp <- enrich.data %>%
      dplyr::filter(cluster == unique(enrich.data$cluster)[x])

    # =============================================
    # enrich

    # entriz id transformation
    id.trans <- ifelse(type == "ownSet",FALSE,TRUE)

    if(id.trans == TRUE){
      gene.ent <- clusterProfiler::bitr(tmp$gene,
                                        fromType = fromType,
                                        toType = toType,
                                        OrgDb = OrgDb)

      tartget.gene = unlist(gene.ent[,c(toType)])
    }else{
      tartget.gene = tmp$gene
    }

    # GO enrich
    if(type %in% c("BP","MF","CC")){
      set.seed(seed)
      ego <- clusterProfiler::enrichGO(gene          = tartget.gene,
                                       keyType       = toType,
                                       OrgDb         = OrgDb,
                                       ont           = type,
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 1,
                                       qvalueCutoff  = 1,
                                       readable      = readable)
    }else if(type == "ownSet"){
      set.seed(seed)
      ego <- clusterProfiler::enricher(gene = tartget.gene,
                                       TERM2GENE     = TERM2GENE,
                                       TERM2NAME     = TERM2NAME,
                                       pvalueCutoff  = 1,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff  = 1)
    }else{
      set.seed(seed)
      ego <- clusterProfiler::enrichKEGG(gene          = tartget.gene,
                                         keyType       = "kegg",
                                         organism      = organism,
                                         universe      = NULL,
                                         pvalueCutoff  = 1,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff  = 1)

      # transform gene id
      if(readable == TRUE){
        ego <- clusterProfiler::setReadable(ego,OrgDb = OrgDb,keyType = toType)
      }else{
        ego <- ego
      }
    }

    # to data.frame
    enrich_res <- data.frame(ego)
    if(nrow(enrich_res) > 0){
      df <- data.frame(ego) %>%
        dplyr::filter(pvalue < pvalueCutoff) %>%
        dplyr::mutate(group = paste("C",unique(enrich.data$cluster)[x],sep = '')) %>%
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
        if(add.gene == TRUE){
          df2 <- df1 %>% dplyr::select(group,Description,pvalue,ratio,geneID)
        }else{
          df2 <- df1 %>% dplyr::select(group,Description,pvalue,ratio)
        }
      }else{
        df2 <- df
      }

      # remove no pass threshold terms
      df2 <- df2 %>% stats::na.omit()

      # results
      return(df2)
    }else{
      return(NULL)
    }
  })
}

#' This is a test data for this package
#' test data describtion
#'
#' @name net
#' @docType data
#' @author Junjun Lao
"net"
