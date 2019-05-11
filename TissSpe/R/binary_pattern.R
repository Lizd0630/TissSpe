#' Ranking psi according to equal-width-intervals
#'
#' Function to cut data.frame range into n equal-width interval points, maximal
#' point and minimal point, then rank them for their value. This function is
#' design for processing data, in which detecting differantial expression rely
#' on diffenrence, like PSI(delta psi).
#'
#' @param df data.frame.
#' @param n number of medium intervals to generated (all \code{n}+2). Default
#' 10.
#' @param min psi under \code{min} will be set to 0.
#' @param max psi over \code{max} will be ranked highest.
#' @return data.frame contain ranks from \code{df}.
psi_seq_rank <- function(df,
                         n = 10,
                         min = 0,
                         max = 100) {
  if (min < 0 | max > 100) {
    stop("Please set min and max psi to 0-100!")
  } else {
    bks <- quantile(0:100, probs = seq(0, 1, 1/n))
  }

  if (min > 0) {
    bks[1] <- min
  } else {
    bks[1] <- bks[1] + 1e-6
  }

  if (max < 100) {
    bks[n+1] <- max
  } else {
    bks[n+1] <- bks[n+1] - 1e-6
  }

  bks <- c(0, bks, 100)

  df_rank <- t(apply(df, 1, function(i) {as.numeric(as.vector(cut(i, breaks = bks, labels = c(0:(n+1)), include.lowest = T)))}))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}



#' Ranking gene expression according to equal-width-intervals
#'
#' Function to cut data.frame range into n equal-width-interval points, maximal
#' point and minimal point, then rank them for their value. This function is
#' design for processed data, like log2 transformed gene expression
#' value(RPKM/FPKM/TPM).
#'
#' @param df data.frame.
#' @param n number of medium intervals to generated (all \code{n}+2). Default 12.
#' @param min psi under \code{min} will be set to 0. Default 0.
#' @param step interval width. Default 1, LFC: log2(2) = 1.
#' @return data.frame contain ranks from \code{df}.
expr_seq_rank <- function(df,
                          n = 12,
                          min = 0,
                          step = 1) {
  max = min + step * n
  df_min <- min(df, na.rm = TRUE)
  df_max <- max(df, na.rm = TRUE)
  if (min < df_min | max > df_max) {
    stop(paste0("Range of expression is ", df_min, "-", df_max, "! Please check parameters!"))
  } else {
    bks <- seq(min, max, step)
    # bks <- quantile(min:max, probs = seq(0, 1, 1/(n-1)))
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

  df_rank <- t(apply(df, 1, function(i) {as.numeric(as.vector(cut(i, breaks = bks, labels = c(0:(n+1)), include.lowest = T)))}))
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
#' @param min psi under \code{min} will be set to 0. Default 0.
#' @param max interval width. Default 12, 2^12 = 4096.
#' @return data.frame contain ranks from \code{df}.
expr_quant_rank <- function(df,
                            n = 10,
                            min = 0,
                            max = 12) {
  vect <- as.vector(t(df))
  df_min <- min(vect, na.rm = TRUE)
  df_max <- max(vect, na.rm = TRUE)

  if (min < df_min | max > df_max) {
    stop(paste0("Please set min and max psi to 0-", df_max, "!"))
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

  df_rank <- t(apply(df, 1, function(i) {as.numeric(as.vector(cut(i, breaks = bks, labels = c(0:(n+1)), include.lowest=T)))}))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}



#' binary index
#'
#' Function to generate binary index(0/1) of ranked data. Define gaps as the
#' diferrence of sorted vectors(ranks, low to high), and ranks over the maximal
#' gap set to 1, otherwise set to 0, always select the first maximal gap. If
#' specify \code{mingap}, then ranks over the \code{mingap} set to 1, othewise
#' set to 0. Ranks value from function:\code{psi_seq_rank}, \code{expr_seq_rank},
#' \code{expr_quant_rank}. Ib is the sum of binary pattern if vector has mingap,
#' otherwise Ib = 0.
#'
#' @param x integer vector, Ranks.
#' @param mingap minimal gap to generate binary pattern.
#' @return binary index with "DE"(differential exprresion) or "UC"(unclear).
bin_index <- function(x,
                      mingap = 3) {
  x_diff <- diff(sort(x))
  if (any(x_diff >= 3)) {
    max_index <- which(x_diff >= 3)[1]
    Ib <- length(x) - max_index
    #cutoff <- x[max_index + 1]
    #x[x < cutoff] <- 0
    #x[x >= cutoff] <- 1
    return(c(Ib, "DE"))
  } else {
    return(c(NA, "UC"))
  }
}



#' binary pattern
#'
#' Function to generate binary pattern(0/1) of ranked data. Define gaps as the
#' diferrence of sorted vectors(ranks, low to high), and ranks over the
#' maximal gap set to 1, otherwise set to 0, always select the first maximal
#' gap. If specify \code{mingap}, then ranks over the \code{mingap} set to 1,
#' othewise set to 0. Ranks value from function:\code{psi_seq_rank},
#' \code{expr_seq_rank}, \code{expr_quant_rank}.
#'
#' @param x integer vector, Ranks.
#' @param mingap minimal gap to generate binary pattern.
#' @return binary pattern, which can be classified into "DE"(differential
#' exprresion, with mingap) or "UC"(unclear, ) or "HK"(house keeping, all in
#' the same rank).
bin_pattern <- function(x,
                      mingap = 3) {
  x_diff <- diff(sort(x))
  if (any(x_diff >= 3)) {
    max_index <- which(x_diff >= 3)[1]
    Ib <- length(x) - max_index
    cutoff <- x[max_index + 1]
    x[x < cutoff] <- 0
    x[x >= cutoff] <- 1
    return(c(x, "DE"))
  } else if (length(unique(x)) == 1) {
    return(c(rep(1, length(x)), "HK"))
  } else {
    return(c(rep(0, length(x)), "UC"))
  }
}
