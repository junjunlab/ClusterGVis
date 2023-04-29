#' @name filter.std
#' @title using filter.std to filter low expression genes
#'
#' @param eset expression matrix, default NULL.
#' @param min.std min stand error, default 0.
#' @param visu whether plot, default FALSE.
#'
#' @return matrix.
filter.std <- function(eset = NULL,
                       min.std = 0,
                       visu = FALSE){
  #index <- logical(dim(exprs(eset))[1])
  tmp <- logical(dim(eset)[1])
  if (is.numeric(min.std)){
    data <- eset
    for (i in 1:length(tmp)){
      tmp[i]   <- stats::sd(data[i,],na.rm=TRUE)
      #   index[i]  <- ( tmp[i] > min.std)

    }
    index <- tmp > min.std
    index[is.na(index)] <- TRUE
    cat(paste(sum(!index),"genes excluded.\n"))
  }

  if (visu){
    plot(sort(tmp),xlab = "Ordered genes",ylab = "Sd")
  }
  eset[index,]
}
