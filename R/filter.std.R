#' @name filter.std modified by Mfuzz filter.std
#' @title using filter.std to filter low expression genes
#'
#' @param eset expression matrix, default NULL.
#' @param min.std min stand error, default 0.
#' @param visu whether plot, default FALSE.
#'
#' @return matrix.
filter.std <- function (eset, min.std,visu=TRUE){
  #index <- logical(dim(exprs(eset))[1])
  if(class(eset) %in% c("data.frame", "matrix")){
    tmp <- logical(dim(eset)[1])
  }else{
    tmp <- logical(dim(Biobase::exprs(eset))[1])
  }

  if (is.numeric(min.std)){
    if(class(eset) %in% c("data.frame", "matrix")){
      data <- eset
    }else{
      data <- Biobase::exprs(eset)
    }

    for (i in 1:length(tmp)){
      tmp[i]  <- sd(data[i,],na.rm=TRUE)
      #   index[i]  <- ( tmp[i] > min.std)
    }
    index <- tmp > min.std
    index[is.na(index)] <- TRUE
    cat(paste(sum(!index),"genes excluded.\n"))
  }


  if (visu){
    plot(sort(tmp),xlab="Ordered genes",ylab="Sd")
  }
  eset[index,]
}
