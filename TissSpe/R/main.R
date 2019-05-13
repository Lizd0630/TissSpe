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
  len <- 1:length(colnames(df))
  df$Tau <- apply(df[, len], 1, ts_Tau)
  df$Gini <- apply(df[, len], 1, ts_Gini)
  df$Tsi <- apply(df[, len], 1, ts_Tsi)
  df$Counts <- apply(df[, len], 1, function(x) {x <- ts_Counts(x, cutoff)})
  df$Hg <- apply(df[, len], 1, ts_Hg)
  df$Zscore <- ts_Z(df[, len])
  df$Spm <- apply(df[, len], 1, ts_Spm)
  df$Ee <- ts_Ee(df[, len])
  df$Pem <- ts_Pem(df[, len])

  df$Mean <- apply(df[, len], 1, expr_mean)
  df$Max <- apply(df[, len], 1, expr_max)
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
#' Function to find tissue-specific AS-events for a given data.frame.
#'
#' Function to detect tissue-specific AS events, and return a list with 2
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
#' @param cutoff Numeric. Values under \code{cutoff} will set to 0(unexpressed)
#' in Specificity method "Counts".
#' @param mingap Integer. Minimal gap of generating binary pattern, Please
#' refer the paper. Default 3.
#' @importFrom stats na.omit
#' @return List of two data.frame. one of them contains psi values with their
#' specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm"(named "raw"), the other contains binary pattern values and
#' index binary "Ib"(named "bin").
#' @export
#' @examples
#' ts_psi(demo_psi,
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

  ## binary type
  df_list <- list(raw = df, bin = psi_seq_rank(df = df, n = n, min = min, max = max))

  ## calculate tissue specificity
  df_list[[1]] <- ts_index(df_list[[1]], cutoff = cutoff)
  df_list[[2]] <- ts_bin(df_list[[2]], mingap = mingap)
  return(df_list)
}



#############################
### For gene expression
#############################
#' Calculate specificity for gene expression
#'
#' Function to find tissue-specific gene expression for a given data.frame. For
#' equal-density-intervals (\code{binary} = "quant"), parameters: \code{df,
#' binary, n, min, max, tissues, identifier} must be specified. However, for
#' equal-width-intervals (\code{binary} = "seq"), parameters: \code{df, binary,
#' n, min, step, tissues, identifier} must be specified.
#'
#' Function to detect tissue specific gene, and return a list with 2
#' data.frames, which contain raw values and binary pattern values and their
#' specificty values of 10 methods: "Tau", "Gini", "Tsi", "Counts", "Ee", "Pem",
#' "Hg", "Z", "Spm", "Ib". Rows with NA will be dropped. Of the binary pattern,
#' all values betwwen \code{min} and \code{max} in the data.frame will be
#' graded into \code{n} equal-width-intervals or \code{n} equal-density-intervals
#' and then assign the rank 0(unexpressed) to \code{n+1}(highest), respectively.
#'
#' @param df data.frame, which contain psi vaules. One column is the names of
#' symbols, like gene id, etc.
#' @param binary "seq" or "quant". Binary-intervals method, "seq" refer to
#' equal-width-intervals and "quant" refer to equal-density-intervals.
#' @param n Integer. \code{n+2} "seq" or "quant" intervals be generated of all
#' expression values. Default 10.
#' @param min Numeric. The values under \code{min} will be graded into rank
#' 0(unexpressed). Be careful when used with \code{trans}. Defaut 0.
#' @param max Numeric. The values over \code{max} will be graded into rank
#' \code{n+1}(highest). Default 16.
#' @param step Numeric. Width of intervals in \code{binary} "seq" method. Be
#' careful when used with \code{trans}.
#' @param trans Charactor. one of "log2", "log2_QN", "QN", "none". "QN" means
#' \code{quantile.normolize}.
#' @param tissues Vector of charactors, at leat length of 2. Analysed Tissues'
#' unique identifiers, and must exactly keep away from "Tau", "Gini", "Tsi",
#' "Counts", "Ee", "Pem", "Hg", "Z", "Spm", "Ib", "Type", "Mean" and "Max".
#' @param identifier Charactor. Length of 1. The colname of unique identifiers
#' for records, like "gene_id", "AS_events", etc.
#' @param na.del Logical. Whether NAs will be dropped. Default TRUE.
#' @param cutoff Numeric. Values under \code{cutoff} will set to 0(unexpressed)
#' in Specificity method "Counts".
#' @param mingap Integer. Minimal gap of generating binary pattern. Default 2.
#' @return List of two data.frame. one of them contains expression values with
#' their specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee",
#' "Pem", "Hg", "Z", "Spm"(named "raw"), the other contains binary pattern values
#' and index binary "Ib"(named "bin").
#' @export
#' @examples
#' ts_psi(demo_tpm,
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

  ## calculate mean of replicates
  df <- rep_mean(df = df, tissues = tissues)

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
  } else if (trans == "none") {
    df  <- df
  } else {
    stop("Value of trans error!")
  }

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
  df_list[[2]] <- ts_bin(df_list[[2]], mingap = mingap)
  return(df_list)
}


