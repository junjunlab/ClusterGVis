#' @name filter.std modified by Mfuzz filter.std
#' @title using filter.std to filter low expression genes
#'
#' @param eset expression matrix, default NULL.
#' @param min.std min stand error, default 0.
#' @param visu whether plot, default FALSE.
#' @param verbose show filter information.
#'
#' @return matrix.
filter.std <- function(eset, min.std, visu = TRUE, verbose = TRUE) {
  # index <- logical(dim(exprs(eset))[1])

  if ("matrix" %in% class(eset) | "data.frame" %in% class(eset)) {
    tmp <- logical(dim(eset)[1])
  } else {
    tmp <- logical(dim(Biobase::exprs(eset))[1])
  }

  if (is.numeric(min.std)) {
    if ("data.frame" %in% class(eset) | "matrix" %in% class(eset)) {
      data <- eset
    } else {
      data <- Biobase::exprs(eset)
    }

    for (i in seq_len(length(tmp))) {
      tmp[i] <- sd(data[i, ], na.rm = TRUE)
      #   index[i]  <- ( tmp[i] > min.std)
    }
    index <- tmp > min.std
    index[is.na(index)] <- TRUE
    if (verbose) {
      methods::show(paste(sum(!index), "genes excluded.\n"))
    }
  }


  if (visu) {
    plot(sort(tmp), xlab = "Ordered genes", ylab = "Sd")
  }
  eset[index, ]
}
