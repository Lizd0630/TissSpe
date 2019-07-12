############################
# binary pattern and index #
############################
#' Calculate binary pattern for a given ranks data.frame
#'
#' Function to calculate binary pattern(0/1) for a given ranks data.frame and
#' return a data.frame with ranks with binary pattern and phenotype("DE" or "UC").
#'
#' @param df data.frame of ranks.
#' @param mingap integer. Minimal gap to generate binary pattern. Default 3.
#' @return data.frame with ranks with binary pattern and phenotype("DE" or "UC").
ts_pattern <- function(df,
                       mingap = 5) {
  res <- as.data.frame(t(apply(df, 1, function(x) {
                               bin_pattern(x, mingap = mingap)
                               })))
  colnames(res) <- c(colnames(df), "Ib", "Type")
  rownames(res) <- rownames(df)
  for (i in 1:(ncol(res) - 1)) {
    res[, i] <- as.integer(as.vector(res[, i]))
  }
  return(res)
}



#' binary pattern
#'
#' Function to generate binary pattern(0/1) of ranked data. Define gaps as the
#' diferrence of sorted vectors(ranks, low to high), and ranks over the
#' maximal gap set to 1, otherwise set to 0, always select the first maximal
#' gap.
#' If specify \code{mingap}, then ranks over the maximal gap set to 1,
#' othewise set to 0. And assign it as "DE", differential expression.
#' If all values in \code{x} are identical, then all binary set to 1. And assign
#' it as "HK", house-keeping.
#' If there is no \code{mingap}, all binary set to 0. And assign it as "UC", unclear.
#' Ranks value from function:\code{psi_seq_rank}, \code{expr_fold_rank},
#' \code{expr_quant_rank}, \code{expr_bks_rank}.
#'
#' @param x integer vector, Ranks.
#' @param mingap minimal gap to generate binary pattern.
#' @return binary pattern, which can be classified into "DE"(Differential
#' exprresion, with mingap) or "UC"(Unclear, without mingap and have different
#' rank) or "HK"(House keeping, all in the same rank).
bin_pattern <- function(x,
                        mingap = 5) {
  if (!is.numeric(x)) {
    stop("x must be numeric!")
  }

  x_sort <- sort(x)
  x_diff <- diff(x_sort)
  if (max(x_diff) == 0) {
    return(c(rep(1, length(x)), length(x), "HK"))
  } else if (any(x_diff >= mingap)) {
    max_index <- which.max(x_diff)
    Ib <- length(x) - max_index
    cutoff <- x_sort[max_index + 1]
    x[x < cutoff] <- 0
    x[x >= cutoff] <- 1
    return(c(x, Ib, "mingap"))
  } else {
    return(c(rep(NA, length(x)), 0, "UC"))
  }
}



#' Ranking psi according to equal-width-intervals
#'
#' Function to cut data.frame range into n equal-width interval points, maximal
#' point and minimal point, then rank them for their value. This function is
#' design for processing data, in which detecting differantial expression rely
#' on diffenrence, like PSI(delta psi).
#'
#' @param df data.frame.
#' @param n number of medium intervals to generated (all \code{n}+2). Default
#' 20.
#' @param min psi under \code{min} will be set to 0.
#' @param max psi over \code{max} will be ranked highest.
#' @return data.frame contain ranks from \code{df}.
psi_seq_rank <- function(df,
                         n = 20,
                         min = 0,
                         max = 100) {
  if (!all(apply(df, 2, is.numeric))) {
    stop("df must be numeric!")
  }

  if (is.numeric(n)) {
    n <- as.integer(n)
  } else {
    stop("n should be integer!")
  }

  if (!is.numeric(min) | min < 0 | min > 100) {
    stop("min should be numeric and range 0-100!")
  }

  if (!is.numeric(max) | max > 100 | max < 0) {
    stop("max should be numeric and range 0-100!")
  }

  bks <- quantile(0:100, probs = seq(0, 1, 1/n))

  if (min > 0 & min < bks[2]) {
    bks[1] <- min
  } else if (min > 0 & min > bks[2]) {
    stop(paste("min should lower than", bks[2]))
  } else {
    bks[1] <- bks[1] + 1e-10
  }

  if (max < 100 & max > bks[n]) {
    bks[n+1] <- max
  } else if (max < 100 & max < bks[n]) {
    stop(paste("max should greater than", bks[n]))
  } else {
    bks[n+1] <- bks[n+1] - 1e-10
  }

  bks <- c(0, bks, 100)

  df_rank <- t(apply(df, 1, function(i) {
    as.numeric(as.vector(cut(i,
                             breaks = bks,
                             labels = c(0:(n+1)),
                             include.lowest = T,
                             rigth = FALSE)))
    }))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}



#' Ranking gene expression according to fold-change-intervals
#'
#' Function to cut data.frame range into \code{2n+2} fold-change-interval
#' points (1,2,3,4,6,8,12,16,24,32,...) and , maximal point and minimal
#' point, then rank them for their value. This function is design for
#' processed data, like log2 transformed gene expression
#' (0,1,1+log(3/2),2,2+log2(3/2)...)value(RPKM/FPKM/TPM).
#'
#' @param df data.frame.
#' @param n number * 2 of medium intervals to generated (all \code{2n}+2).
#' Default 12. Values greater than n set to highest rank.
#' @param min expression under \code{min} will be set to 0. Default 0.5.
#' @param step interval width. Default 1, log2 fold change: log2(2) = 1.
#' @param medium the medium interval of neighbour 2 steps. Default log2(3/2).
#' @return data.frame contain ranks from \code{df}.
expr_fold_rank <- function(df,
                          n = 12,
                          min = 0.5,
                          medium = log2(3/2),
                          step = 1) {
  if (!all(apply(df, 2, is.numeric))) {
    stop("df must be numeric!")
  }

  if (is.numeric(n)) {
    n <- as.integer(n)
  } else {
    stop("n should be integer!")
  }

  if (!is.numeric(min) | min < 0) {
    stop("min should be numeric and greater than 0!!")
  }

  if (!is.numeric(step) | step <= 0) {
    stop("step should be numeric and greater than 0!")
  }

  df_min <- min(df, na.rm = TRUE)
  df_max <- max(df, na.rm = TRUE)
  if ((step * n) > df_max) {
    stop(paste("The maximum n you can set is", as.integer(df_max/step), "!"))
  } else {
    message(paste("The maximum n you can set is", as.integer(df_max/step), "!"))
  }

  bks <- sort(c(seq(0, step * n, step), (seq(1, step * n, step) + medium)))
  m <- length(bks)

  if (min == 0 | min >= bks[2]) {
    bks[1] <- 1e-10
  } else if (min < bks[2]) {
    bks[1] <- min
  }

  if (max(bks) > df_max) {
    bks <- bks[-m]
  } else if ((max(bks) == df_max)) {
    bks[m] <- bks[m] - 1e-10
  }

  bks <- c(0, bks, df_max)


  df_rank <- t(apply(df, 1, function(i) {
    as.numeric(as.vector(cut(i,
                             breaks = bks,
                             labels = c(0:(length(bks)-2)),
                             include.lowest = T,
                             right = FALSE)))
    }))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}



#' Ranking gene expression according to equal-density-intervals(quantiles)
#'
#' Function to cut data.frame range into n equal-density-interval(quantiles)
#' points, maximal point and minimal point, then rank them for their value.
#' This function is design for processed data, like log2 transformed gene
#' expression value(RPKM/FPKM/TPM).
#'
#' @param df data.frame.
#' @param n number of medium intervals to generated (all \code{n}+2). Default
#' 10.
#' @param min expression under \code{min} will be set to 0. Default 0.
#' @param max interval width. Default 12, 2^12 = 4096.
#' @return data.frame contain ranks from \code{df}.
expr_quant_rank <- function(df,
                            n = 10,
                            min = 0,
                            max = 12) {
  if (!all(apply(df, 2, is.numeric))) {
    stop("df must be numeric!")
  }

  vect <- as.vector(t(df))
  df_min <- min(vect, na.rm = TRUE)
  df_max <- max(vect, na.rm = TRUE)

  if (min < df_min | max > df_max) {
    stop(paste0("Please set min and max to 0-", df_max, "!"))
  } else {
    vect <- vect[(vect >= min & vect <= max)]
    bks <- as.vector(quantile(vect, probs = seq(0, 1, 1/n)))
  }

  if (min > df_min) {
    bks[1] <- min
  } else {
    bks[1] <- bks[1] + 1e-6
  }

  if (max < df_max) {
    bks[n+1] <- max
  } else {
    bks[n+1] <- bks[n+1] - 1e-6
  }

  bks <- c(df_min, bks, df_max)

  df_rank <- t(apply(df, 1, function(i) {
    as.numeric(as.vector(cut(i,
                             breaks = bks,
                             labels = c(0:(n+1)),
                             include.lowest=T,
                             right = F)))
    }))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}



#' Ranking gene expression according to given breaks
#'
#' Function to cut data.frame range intointervals according to given breaks
#' points, maximal point and minimal point, then rank them for their value.
#' This function is design for raw expression value(RPKM/FPKM/TPM).
#'
#' @param df data.frame.
#' @param bks vector, numeric. Intervals to cut according to \code{bks}.
#' Finally generate \code{length(bks)-1} ranks.
#' @return data.frame contain ranks from \code{df}.
expr_bks_rank <- function(df,
                          bks) {
  if (!all(apply(df, 2, is.numeric))) {
    stop("df must be numeric!")
  }

  if (!is.numeric(bks)) {
    stop("bks should be numeric vectors!")
  }

  bks <- sort(bks)
  min <- min(bks)
  max <- max(bks)
  n <- length(bks)

  vect <- as.vector(t(df))
  df_min <- min(vect, na.rm = TRUE)
  df_max <- max(vect, na.rm = TRUE)


  if (min > df_min | max < df_max) {
    stop(paste0("bks should range ", df_min, "-", df_max, "!"))
  }

  df_rank <- t(apply(df, 1, function(i) {
    as.numeric(as.vector(cut(i,
                             breaks = bks,
                             labels = c(1:(n-1)),
                             include.lowest = T,
                             right = F)))
    }))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}


