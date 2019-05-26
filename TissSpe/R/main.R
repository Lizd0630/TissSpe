#########################
###   access tissues  ###
#########################
#' check whether exsiting substring and return messages.
#' @param vect vector th access.
chk_sub <- function(vect) {
  mes <- c()
  for (i in 1:length(vect)) {
    rest <- vect[-i]
    hits <- grepl(paste0("^", vect[i]), rest)
    if (any(hits)) {
      mes <- c(mes, paste("tissues", vect[i], "is the substring of ",
                          paste(rest[hits], collapse = ",")))
    }
  }
  return(mes)
}



#########################
### For AS events psi ###
#########################
#' Calculate specificity of psi
#'
#' Function to find tissue-specific AS-events for a given data.frame.
#'
#' Function to detect tissue-specific AS events, and return a list with 3
#' data.frames, which contain psi values and binary pattern values and their
#' specificty values of 10 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm", "Ib". Rows with NA will be dropped. Of the binary pattern,
#' all values betwwen min and max in the data.frame will be graded into n
#' equal-width-intervals and then assign the rank 0(unexpressed) to \code{n+1}.
#'
#' @param df data.frame, which contain psi vaules (0-100). One column is the names of
#' symbols, like AS events id, etc.
#' @param n Integer. \code{n+2} equal-width-intervals be generated of all psi
#' values. Default 10.
#' @param min Numeric. The values under \code{min} will be graded into rank
#' 0(unexpressed). Defaut 0.
#' @param max Numeric. The values over \code{max} will be graded into rank
#' \code{n+1}(highest). Default 100.
#' @param tissues Vector of charactors, at leat length of 2. Analysed
#' Tissues' unique identifier, and must keep away from "Tau", "Gini", "Tsi",
#' "Counts", "Ee", "Pem", "Hg", "Z", "Spm", "Ib", "Type", "Mean" and "Max",
#' exactly.
#' @param identifier Charactor, length of 1. The colname of unique identifier
#' for row symbols, like "gene_id", "AS_events", etc.
#' @param na.del Logical. Whether NAs will be dropped. Default TRUE.
#' @param mingap Integer. Minimal gap of generating binary pattern, Please
#' refer the paper. Default 5.
#' @importFrom stats na.omit
#' @return List of 3 data.frame. one of them contains psi values with their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm"(named "raw"), the second contains rank values and binary
#' index(named "rank"), the third contains binary pattern values and
#' index binary "Ib"(named "bin").
#' @export
#' @examples
#' res1 <- ts_psi(demo_psi, n = 20, min = 0.5,
#'                tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                            "sample_E", "sample_F", "sample_G", "sample_H",
#'                            "sample_I", "sample_J", "sample_K", "sample_L",
#'                            "sample_M", "sample_N", "sample_O", "sample_P",
#'                            "sample_Q"),
#'                            identifier = "AS_events")
#' res2 <- ts_psi(demo_psi, n = 10, min = 1,
#'                tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                            "sample_E", "sample_F", "sample_G", "sample_H",
#'                            "sample_I"),
#'                            identifier = "AS_events")
ts_psi <- function(df,
                  n = 20,
                  min = 0,
                  max = 100,
                  tissues,
                  identifier,
                  na.del = TRUE,
                  mingap = 3) {
  ## check tissue names
  mes <- chk_sub(tissues)
  if (!is.null(mes)) {
    stop(mes[1])
  }

  if (!is.numeric(min) | min < 0 | min > 100) {
    stop("min should be numeric, 0-100!")
  }

  if (!is.numeric(max) | max < 0 | max > 100) {
    stop("max should be numeric, 0-100!")
  }

  if (!is.numeric(n) | n <= 0) {
    stop("n should be integer and greater than 0!")
  } else {
    n <- as.integer(n)
  }

  if (!is.numeric(mingap) | mingap <= 0) {
    stop("mingap should be integer and greater than 0!")
  } else {
    mingap <- as.integer(mingap)
  }

  ## format data.frame
  if (is.vector(tissues) & length(tissues) >= 2) {
    df <- fmt_df(df = df, tissues = tissues, identifier = identifier)
    if (max(df) <= 1) {
      stop("Vaules of psi should convert into 0-100!")
    }
  } else {
    stop("tissues must be a vector with length at least 2!")
  }

  ## calculate mean of replicates
  df <- rep_mean(df = df, tissues = tissues)

  ## wheather remove NA
  if (na.del == TRUE) {
    df <- na.omit(df)
  }

  df[df < min] <- 0

  ## binary type
  df_list <- list(raw = df, rank = psi_seq_rank(df = df, n = n, min = min, max = max))

  ## calculate tissue specificity and binary index and binary pattern
  df_list$raw <- ts_index(df_list$raw, cutoff = min)
  df_list$bin <- ts_pattern(df_list$rank, mingap = mingap)
  df_list$rank$Ib <- df_list$bin$Ib
  df_list$rank$Type <- df_list$bin$Type
  return(df_list)
}



#############################
### For gene expression
#############################
#' Calculate specificity for gene expression
#'
#' Function to find tissue-specific gene expression for a given data.frame.
#' For equal-density-intervals (\code{binary = "quant"}), parameters: \code{df,
#' binary, n, min, max, tissues, identifier} must be specified, and the used
#' values are in the range \code{(min, max)}, values lower than min class as 0,
#' and higher than max class as highest rank.
#' However, for fold-intervals (\code{binary = "fold"}), parameters: \code{df,
#' binary, n, min, tissues, identifier} must be specified, and the used values
#' are in the range \code{(min, n)}.
#' Finally, for given breaks (\code{binary = "bks"}), parameters: \code{df,
#' binary, bks, tissues, identifier} must be specified.
#'
#' Function to detect tissue specific gene, and return a list with 2
#' data.frames, which contain raw values and binary pattern values and their
#' specificty values of 10 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm", "Ib". Rows with NA will be dropped. Of the binary pattern,
#' all values betwwen \code{min} and \code{max} in the data.frame will be
#' graded into fold-intervals or \code{n} equal-density-intervals
#' and then assign the rank 0(unexpressed) to \code{n+1}(highest), respectively.
#'
#' @param df data.frame, which contain expression vaules. One column is the names
#' of symbols, like gene_id, etc.
#' @param binary "fold", "quant" or "bks". Binary-intervals method, "fold" refer to
#' fold-change-intervals and "quant" refer to equal-density-intervals. "bks"
#' refer to given breaks points interval.
#' @param n Integer. \code{n+2} "fold" or "quant" intervals be generated of all
#' expression values. Default 10.
#' @param min Numeric. The values under \code{min} will be graded into rank
#' 0(unexpressed). Be careful when used with \code{trans}. Defaut 0.
#' @param max Numeric. The values over \code{max} will be graded into rank
#' \code{n+1}(highest). Default 16.
#' @param bks Vector, numercic. Design for method \code{binary="bks"}. Ranks will
#' assigned according to \code{bks}.
#' @param trans Charactor. one of "log2", "log2_QN", "QN", "none". "QN" means
#' \code{quantile.normolize}.
#' @param tissues Vector of charactors, at leat length of 2. Analysed Tissues'
#' unique identifiers, and must exactly keep away from "Tau", "Gini", "Tsi",
#' "Counts", "Ee", "Pem", "Hg", "Z", "Spm", "Ib", "Type", "Mean" and "Max".
#' @param identifier Charactor. Length of 1. The colname of unique identifiers
#' for records, like "gene_id", "AS_events", etc.
#' @param na.del Logical. Whether NAs will be dropped. Default TRUE.
#' @param mingap Integer. Minimal gap of generating binary pattern. Default 3.
#' @return List of 3 data.frame. one of them contains psi values with their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm"(named "raw"), the second contains rank values and binary
#' index(named "rank"), the third contains binary pattern values and
#' index binary "Ib"(named "bin").
#' @export
#' @examples
#' res1 <- ts_expr(demo_tpm, binary = "fold", n = 11, min = 0.5,
#'                 tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                             "sample_E", "sample_F", "sample_G", "sample_H",
#'                             "sample_I", "sample_J", "sample_K", "sample_L",
#'                             "sample_M", "sample_N", "sample_O", "sample_P"),
#'                             identifier = "gene_id")
#' res2 <- ts_expr(demo_tpm, binary = "quant", min = 0.5, max = 11,
#'                 tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                             "sample_E", "sample_F", "sample_G", "sample_H",
#'                             "sample_I", "sample_J", "sample_K", "sample_L",
#'                             "sample_M", "sample_N", "sample_O", "sample_P"),
#'                             identifier = "gene_id")
#' res3 <- ts_expr(demo_tpm, binary = "bks", min = 0.5, max = 11,
#'                 bks = seq(0,18,1),
#'                 tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                             "sample_E", "sample_F", "sample_G", "sample_H",
#'                             "sample_I", "sample_J", "sample_K", "sample_L",
#'                             "sample_M", "sample_N", "sample_O", "sample_P"),
#'                             identifier = "gene_id")
ts_expr <- function(df,
                    binary = "fold",
                    n = 12,
                    min = 0,
                    max = 16,
                    bks,
                    trans = "log2",
                    tissues,
                    identifier,
                    na.del = TRUE,
                    mingap = 3) {
  ## check tissue names
  mes <- chk_sub(tissues)
  if (!is.null(mes)) {
    stop(mes[1])
  }

  if (!is.numeric(min) | min < 0 | min > 100) {
    stop("min should be numeric, 0-100!")
  }

  if (!is.numeric(max) | max < 0 | max > 100) {
    stop("max should be numeric, 0-100!")
  }

  if (!is.numeric(n) | n <= 0) {
    stop("n should be integer and greater than 0!")
  } else {
    n <- as.integer(n)
  }

  if (!is.numeric(mingap) | mingap <= 0) {
    stop("mingap should be integer and greater than 0!")
  } else {
    mingap <- as.integer(mingap)
  }

  ## format data.frame
  if (is.vector(tissues) & length(tissues) >= 2) {
    df <- fmt_df(df = df, tissues = tissues, identifier = identifier)
  } else {
    stop("tissues must be a vector with length at least 2!")
  }

  ## calculate mean of replicates
  df <- rep_mean(df = df, tissues = tissues)

  ## wheather remove NA
  if (na.del == TRUE) {
    df <- na.omit(df)
  }

  cutoff <- min
  if (trans == "log2_QN") {
    df <- quant_norm(log2(df + 1))
    cutoff <- log2(cutoff + 1)
    df[df < cutoff] <- 0
  } else if (trans == "log2") {
    df <- log2(df + 1)
    cutoff <- log2(cutoff + 1)
    df[df < cutoff] <- 0
  } else if (trans == "QN") {
    df <- quant_norm(df)
    df[df < cutoff] <- 0
  } else if (trans == "none") {
    df[df < cutoff] <- 0
  } else {
    stop("Value of trans error!")
  }

  ## binary type
  if (binary == "fold") {
    df_list <- list(raw = df, rank = expr_fold_rank(df = df, n = n, min = min))
  } else if (binary == "quant") {
    df_list <- list(raw = df, rank = expr_quant_rank(df = df, n = n, min = min, max = max))
  } else if (binary == "bks") {
    df_list <- list(raw = df, rank = expr_bks_rank(df = df, bks = bks))
  } else {
    stop("binary error!")
  }

  ## calculate tissue specificity and binary index and binary pattern
  df_list$raw <- ts_index(df_list$raw, cutoff = cutoff)
  df_list$bin <- ts_pattern(df_list$rank, mingap = mingap)
  df_list$rank$Ib <- df_list$bin$Ib
  df_list$rank$Type <- df_list$bin$Type
  return(df_list)
}


