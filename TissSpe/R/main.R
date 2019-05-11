##############################
### For caculate specificity
##############################
#' Calculate specificity for a given numeric data.frame
#'
#' Function to calculate specificity for a given numeric data.frame and return
#' a data.frame with original value with specificity index of 9 methhods.
#'
#' @param df data.frame of numeric.
#' @param cutoff Numeric. Values under cutoff will set to 0(unexpressed) in
#' Specificity method Counts.
#' @return data.frame with original value with specificity index of 9 methhods.
ts_index <- function(df,
                     cutoff = 1) {
  df$Tau <- apply(df, 1, ts_Tau)
  df$Gini <- apply(df, 1, ts_Gini)
  df$Tsi <- apply(df, 1, ts_Tsi)
  df$Counts <- apply(df, 1, function(x) {x <- ts_Counts(x, cutoff)})
  df$Hg <- apply(df, 1, ts_Hg)
  df$Zscore <- ts_Z(df)
  df$Spm <- apply(df, 1, ts_Spm)
  df$Ee <- ts_Ee(df)
  df$Pem <- ts_Pem(df)

  df$Mean <- apply(df, 1, expr_mean)
  df$Max <- apply(df, 1, expr_max)
  return(df)
}



#' Calculate binary index for a given ranks data.frame
#'
#' Function to calculate binary index(0/1) for a given ranks data.frame and
#' return a data.frame with ranks with binary index and phenotype("DE" or "UC").
#'
#' @param df data.frame of ranks.
#' @param mingap integer. Minimal gap to generate binary pattern. Default 3.
#' @return data.frame with ranks with binary index and phenotype("DE" or "UC").
ts_bin <- function(df,
                   mingap = 3) {
  res <- t(apply(df, 1, function(x) {bin_index(x, mingap = mingap)}))
  df$Ib <- as.integer(res[, 1])
  df$Type <- res[, 2]
  return(df)
}



##############################
### For AS events psi
##############################
#' Calculate specificity of psi
#'
#' \code{ts_psi} returns the list of original data with specificities of 10
#' methods.
#'
#' Function to detect tissue-specific AS events, and return a list with 2
#' data.frames, which contain raw values and binary pattern values and their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm". Rows with NA will be dropped. Of the binary pattern, all
#' values betwwen min and max in the data.frame will be graded into n
#' equal-width-intervals and assign the rank 0 to (n+1), respectively.
#'
#' @param df data.frame, which contain psi vaules.
#' @param n Integer. Cut the psi values into (n+2) equal-width-intervals.
#' Default 10.
#' @param min Numeric. The values under min will be graded into rank 0.
#' Defaut 0.
#' @param max Numeric. The values over max will be graded into rank (n+1).
#' Default 100.
#' @param tissues Vector of charactors, at leat 2. Analysed Tissues' unique
#' identifier, and must keep away from "Tau", "Gini", "Tsi", "Counts", "Ee",
#' "Pem", "Hg", "Z", "Spm", "Mean" and "Max".
#' @param identifier Charactor. The colname of unique identifier for row
#' symbols, like "gene_id" or "AS_events".
#' @param na.del Logical. Weather NA will be dropped. Default TRUE.
#' @param cutoff Numeric. Values under cutoff will set to 0(unexpressed) in
#' Specificity method Counts.
#' @param mingap integer. Minimal gap to generate binary pattern. Default 3.
#' @importFrom stats na.omit
#' @return Specificity of psi values in \code{df}. A list with 2 data.frames,
#' which contain raw values and binary pattern values and their specificty
#' values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem", "Hg",
#' "Z", "Spm".
#' @export
#' @examples
#' ts_psi(tmp_psi,
#'        tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                    "sample_E", "sample_F", "sample_G", "sample_H",
#'                    "sample_I", "sample_J", "sample_K", "sample_L",
#'                    "sample_M", "sample_N", "sample_O", "sample_P",
#'                    "sample_Q"),
#'                    identifier = "AS_events")
ts_psi <- function(df,
                  n = 10,
                  min = 0,
                  max = 100,
                  tissues,
                  identifier,
                  na.del = TRUE,
                  cutoff = 1,
                  mingap = 3) {
  ## format data.frame
  if (is.vector(tissues) & length(tissues) >= 2) {
    df <- fmt_df(df = df, tissues = tissues, identifier = identifier)
  } else {
    stop("tissues must be a vector with length at least 2!")
  }

  ## wheather remove NA
  if (na.del == TRUE) {
    df <- na.omit(df)
  }

  df <- rep_mean(df = df, tissues = tissues)

  ## binary type
  df_list <- list(raw = df, bin = psi_seq_rank(df = df, n = n, min = min, max = max))
#  if (binary == "seq") {
#    df_list <- list(raw = df, bin = seq_rank(df = df, n = n, min = min, max = max))
#  } else if (binary == "quant") {
#    df_list <- list(raw = df, bin = quant_rank(df = df, n = n, min = min, max = max))
#  } else {
#    stop("binary type error!")
#  }

  ## calculate tissue specificity
  df_list[[1]] <- ts_index(df_list[[1]], cutoff = cutoff)
  #df_list[[2]] <- ts_index(df_list[[2]], cutoff = cutoff)
  df_list[[2]] <- ts_bin(df_list[[2]], mingap = mingap)
  return(df_list)
}






#############################
### For gene expression
#############################
#' Calculate specificity for a given gene expression data.frame
#'
#' Function to detect tissue specific gene, and return a list with 2
#' data.frames, which contain raw values and binary pattern values and their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm".
#' Rows with NA will be dropped. Of the binary pattern, all values betwwen min
#' and max in the data.frame will be graded into n equal-width-intervals or n
#' equal-density-intervals and assign the rank 0 to (n+1), respectively.
#'
#' @param df data.frame, which contain gene expression vaules.
#' @param binary "seq" or "quant". Binary method, "seq" refer to
#' equal-width-intervals and "quant" refer to equal-density-intervals.
#' @param n Integer. Cut the gene expression values into (n+2)
#' equal-density-intervals. Default 10.
#' @param min Numeric. The values under min will be graded into rank 0.
#' Defaut 0. Be careful when used with \code{trans}.
#' @param max Numeric. The values over max will be graded into rank (n+1).
#' Default 16. Be careful when used with \code{trans}.
#' @param step Numeric. Width of intervals in "seq" method. Be careful when
#' used with \code{trans}.
#' @param trans Charactor. "log2", "log2_QN", "QN". "QN" means
#' \code{quantile.normolize}.
#' @param tissues Vector of charactors, at leat 2. Analysed Tissues' unique
#' identifiers, and must keep away from "Tau", "Gini", "Tsi", "Counts", "Ee",
#' "Pem", "Hg", "Z", "Spm", "Mean" and "Max".
#' @param identifier Charactor. The colname of unique identifiers for row
#' symbols, like "gene_id" or "AS_events".
#' @param na.del Logical. Weather NA will be dropped. Default TRUE.
#' @param cutoff Numeric. Values under cutoff will set to 0(unexpressed) in
#' Specificity method Counts.
#' @param mingap integer. Minimal gap to generate binary pattern. Default 2.
#' @return Specificity of gene expression values in \code{df}. A list with 2
#' data.frames, which contain raw values and binary pattern values and their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm".
#' @export
#' @examples
#' ts_psi(tmp_tpm,
#'        tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                    "sample_E", "sample_F", "sample_G", "sample_H",
#'                    "sample_I", "sample_J", "sample_K", "sample_L",
#'                    "sample_M", "sample_N", "sample_O", "sample_P"),
#'                    identifier = "gene_id")
ts_expr <- function(df,
                    binary = "seq",
                    n = 12,
                    min = 0,
                    max = 16,
                    step = 1,
                    trans = "log2_QN",
                    tissues,
                    identifier,
                    na.del = TRUE,
                    cutoff = 1,
                    mingap = 2) {
  ## format data.frame
  if (is.vector(tissues) & length(tissues) >= 2) {
    df <- fmt_df(df = df, tissues = tissues, identifier = identifier)
  } else {
    stop("tissues must be a vector with length at least 2!")
  }

  ## wheather remove NA
  if (na.del == TRUE) {
    df <- na.omit(df)
  }

  if (trans == "log2_QN") {
    df[df < cutoff] <- 1
    df <- quant_norm(log2(df))
    cutoff <- log2(cutoff)
  } else if (trans == "log2") {
    df[df < cutoff] <- 1
    df <- log2(df)
    cutoff <- log2(cutoff)
  } else if (trans == "QN") {
    df[df < cutoff] <- 0
    df <- quant_norm(df)
  } else {
    df  <- df
  }

  ## calculate mean of replicates
  df <- rep_mean(df = df, tissues = tissues)

  ## binary type
  if (binary == "seq") {
    df_list <- list(raw = df, bin = expr_seq_rank(df = df, n = n, min = min, step = step))
  } else if (binary == "quant") {
    df_list <- list(raw = df, bin = expr_quant_rank(df = df, n = n, min = min, max = max))
  } else {
    stop("binary type error!")
  }

  ## calculate tissue specificity
  df_list[[1]] <- ts_index(df_list[[1]], cutoff = cutoff)
  #df_list[[2]] <- ts_index(df_list[[2]], cutoff = cutoff)
  df_list[[2]] <- ts_bin(df_list[[2]], mingap = mingap)
  return(df_list)
}


