############################
# binary pattern and index #
############################
#' Calculate binary index for a given ranks data.frame
#'
#' Function to calculate binary index(0/1) for a given ranks data.frame and
#' return a data.frame with ranks with binary index and phenotype("DE" or "UC").
#'
#' @param df data.frame of ranks.
#' @param mingap integer. Minimal gap to generate binary pattern. Default 3.
#' @return data.frame with ranks with binary index and phenotype("DE" or "UC").
ts_pattern <- function(df,
                  mingap = 3) {

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



#' binary index
#'
#' Function to generate binary index(the sum of \code{binary pattern}) of
#' ranked data. Define gaps as the diferrence of sorted vector(ranks,
#' low to high), and ranks over the maximal gap in the sorted vector set to 1,
#' otherwise set to 0, always select the first maximal gap.
#' If specify \code{mingap}, then ranks over the \code{mingap} set to 1, othewise
#' set to 0. And if there has no \code{mingap} in the vector, binary index set
#' to NA. Ranks value from function:\code{psi_seq_rank}, \code{expr_seq_rank},
#' \code{expr_quant_rank}. Ib is the sum of binary pattern if vector has mingap,
#' otherwise Ib = 0.
#'
#' @param x Integer vector, Ranks.
#' @param mingap Minimal gap to generate binary pattern.
#' @return binary index with "DE"(differential exprresion) or "UC"(unclear).
bin_index <- function(x,
                      mingap = 3) {
  if (!is.numeric(x)) {
    stop("x must be numeric!")
  }

  if (length(x) <= 1) {
    stop("x must be larger than 1!")
  }

  x_diff <- diff(sort(x))
  if (max(x_diff) == 0) {
    return(c(rep(1, length(x)), "HK"))
  } else if (any(x_diff >= mingap)) {
    max_index <- which.max(x_diff)
    Ib <- length(x) - max_index
    return(c(Ib, "mingap"))
  } else {
    return(c(NA, "UC"))
  }
}



#' binary pattern
#'
#' Function to generate binary pattern(0/1) of ranked data. Define gaps as the
#' diferrence of sorted vectors(ranks, low to high), and ranks over the
#' maximal gap set to 1, otherwise set to 0, always select the first maximal
#' gap.
#' If specify \code{mingap}, then ranks over the \code{mingap} set to 1,
#' othewise set to 0.
#' If all values in \code{x} are identical, then all binary set to 1.
#' If there is no \code{mingap}, all binary set to 0.
#' Ranks value from function:\code{psi_seq_rank}, \code{expr_seq_rank},
#' \code{expr_quant_rank}.
#'
#' @param x integer vector, Ranks.
#' @param mingap minimal gap to generate binary pattern.
#' @return binary pattern, which can be classified into "DE"(Differential
#' exprresion, with mingap) or "UC"(Unclear, without mingap and have different
#' rank) or "HK"(House keeping, all in the same rank).
bin_pattern <- function(x,
                        mingap = 3) {
  if (!is.numeric(x)) {
    stop("x must be numeric!")
  }

  x_diff <- diff(sort(x))
  if (max(x_diff) == 0) {
    return(c(rep(1, length(x)), length(x), "HK"))
  } else if (any(x_diff >= mingap)) {
    max_index <- which.max(x_diff)
    Ib <- length(x) - max_index
    cutoff <- x_diff[max_index + 1]
    x[x < cutoff] <- 0
    x[x >= cutoff] <- 1
    return(c(x, sum(x), "mingap"))
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
#' 10.
#' @param min psi under \code{min} will be set to 0.
#' @param max psi over \code{max} will be ranked highest.
#' @return data.frame contain ranks from \code{df}.
psi_seq_rank <- function(df,
                         n = 10,
                         min = 0,
                         max = 100) {
  if (!all(apply(df, 2, is.numeric))) {
    stop("df must be numeric!")
  }

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
  if (!all(apply(df, 2, is.numeric))) {
    stop("df must be numeric!")
  }

  max = min + step * n
  df_min <- min(df, na.rm = TRUE)
  df_max <- max(df, na.rm = TRUE)
  if (min < df_min | max > df_max) {
    stop(paste0("Range of expression is ", df_min, "-", df_max, "! Please check parameters!"))
  } else {
    warning(paste0("Range of expression is ", df_min, "-", df_max, "! Be careful!"))
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
  if (!all(apply(df, 2, is.numeric))) {
    stop("df must be numeric!")
  }

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



