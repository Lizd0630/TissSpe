#' mean of vectors
#'
#' Function require a vector with numeric and then calculate the
#' mean. Mean is calculated taking in account tissues with 0 expression:
#' 2+0+4=2.
#'
#' @param x Vector, numeric.
expr_mean <- function(x) {
  if(!all(is.na(x))) {
    res <- mean(x, na.rm = T)
  } else {
    res <- NA
  }
  return(res)
}



#' maximal value of replicates
#'
#' Function require a vector with numeric in different conditions.
#' Max is calculated taking in account tissues with 0 expression: 2+0+4=4.
#'
#' @param x numeric.
expr_max <- function(x) {
  if(!all(is.na(x))) {
    res <- max(x, na.rm = T)
  } else {
    res <- NA
  }
  return(res)
}



#' format data.frame
#'
#' Generate a data.frame, which only contain exprssion data of \code{tissues}
#' and rownames is value of \code{identifier} (gene_id, events).
#'
#' @param df numeric.
#' @param tissues charactors. Colnames you want in final results.
#' @param identifier which col to rownames.
fmt_df <- function(df,
                   tissues,
                   identifier) {
  cols <- colnames(df)
  cols <- unique(unlist(sapply(tissues, function(x) {cols[regexpr(x, cols) > 0]})))
  df_new <- df[, cols]
  rownames(df_new) <- df[, identifier]
  return(df_new)
}



#' mean expression of replicates
#'
#' Function requires data frame with numeric. \code{rowMeans} between
#' replicates are calculated. \code{tissues} must be the unique word-start
#' identifier to recognize sample replicates. \code{identifier} is the colname
#' of gene names or AS events.
#'
#' @param df data.frame.
#' @param tissues charactors. Unique word-start identifiers to recognize
#' sample replicates.
#
rep_mean <- function(df,
                     tissues) {
  mat <- matrix(NA, ncol = length(tissues), nrow = nrow(df))
  colnames(mat) <- tissues
  rownames(mat) <- rownames(df)
  for (i in tissues) {
    tissue <- (regexpr(paste0("^", i), colnames(df)) > 0)
    if (!any(tissue)) {
      stop(paste0("'", i, "'", " may not in your data! Please check!"))
    }
    mat[, i] <- rowMeans(as.data.frame(df[, tissue]),
                         na.rm = TRUE,
                         dims = 1)
  }
  df <- as.data.frame(mat)
  return(df)
}



#' quantile normalization
#'
#' Quantile-normalize data.frame. All 0 are set to NA, to exclude them from
#' quantile-normalization, then 0 values (the one set to NA) are set back to
#' 0 finally.
#'
#' @param df data.frame.
#' @importFrom preprocessCore normalize.quantiles
quant_norm <- function(df) {
  df[df == 0] <- NA
  mat <- as.matrix(df)
  mat <- normalize.quantiles(mat)
  mat[is.na(mat)] <- 0
  colnames(mat) <- colnames(df)
  rownames(mat) <- rownames(df)
  return(as.data.frame(mat))
}



#' Sort data.frame (decrease)
#'
#' Sort rows of dataframe according row maximal value coordinate.
#'
#' @param dat data.frame, numeric.
sort_dat_de <- function(dat) {
  if (all(apply(dat, 2, is.numeric))) {
    maxindex <- apply(dat, 1, which.max)
    new_dat <- list()
    for (i in sort(maxindex)) {
      a <- dat[maxindex == i, ]
      a <- a[order(-a[, i]),]
      new_dat[[i]] <- a
    }
    new_dat <- do.call(rbind, new_dat)
    return(new_dat)
  } else {
    stop("dat must be numeric!")
  }
}



#' Sort data.frame (increase)
#'
#' Sort rows of dataframe according row minimal value coordinate.
#'
#' @param dat data.frame, numeric.
sort_dat_in <- function(dat) {
  if (all(apply(dat, 2, is.numeric))) {
    minindex <- apply(dat, 1, which.min)
    new_dat <- list()
    for (i in sort(minindex)) {
      a <- dat[minindex == i, ]
      a <- a[order(a[, i]),]
      new_dat[[i]] <- a
    }
    new_dat <- do.call(rbind, new_dat)
    return(new_dat)
  } else {
    stop("dat must be numeric!")
  }
}


