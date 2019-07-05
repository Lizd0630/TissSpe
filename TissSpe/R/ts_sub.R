#' tissue-specific gene subset
#'
#' Function to get a gene subset accroding to specified \code{Ib} and
#' \code{subtissues}. The input of data is from \code{ts_psi} or \code{ts_expr}.
#'
#' @param lst List of 3 data.frame. one of them contains raw values(or psi values
#' with their specificty values of 9 methods: "Tau", "Gini", "Tsi", "Counts", "Ee",
#' "Pem", "Hg", "Z", "Spm"(named "raw"), the second contains rank values and binary
#' index(named "rank"), the third contains binary pattern values and index binary "Ib"
#' (named "bin"), which were generated from \code{ts_psi} or \code{ts_expr}.
#' @param Ib Integer. Binary index, and is the length of \code{subtissues}.
#' @param ts_method Specificity methods. "Tau", "Gini", "TSI", "Counts", "EE",
#' "Hg", "Zscore", "Spm" or "Pem".
#' @param specificity Vector of range, numeric. Region of specificity. Default
#' \code{c(0, 1)}.
#' @param subtissues Character or NULL. Of binary pattern, character: only tissues in
#' \code{subtissues} have value 1, that is specific-overexpression. NULL: all the
#' Ib will be included.
#' @return a list with 5 data.frames. Subset of \code{lst}.
#' @param sample_mean c(-Inf, Inf)
#' @param sample_mean_sd Inf
#' @param sample_max c(-Inf, Inf)
#' @param sample_range 0
#' @param rep_diff_max Inf
#' @param rep_diff_sd Inf
#' @param rep_diff_mean Inf
#' @export
#' @examples
#' res1 <- ts_expr(demo_tpm, n = 11,
#'         tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                    "sample_E", "sample_F", "sample_G", "sample_H",
#'                    "sample_I", "sample_J", "sample_K", "sample_L",
#'                    "sample_M", "sample_N", "sample_O", "sample_P"),
#'                    identifier = "gene_id")
#' ts_sub(res1, Ib = 2, specificity = c(0, 1),
#'        subtissues = c("sample_B", "sample_C"))
#' ts_sub(res1, Ib = 3)
#' res2 <- ts_psi(demo_psi,
#'         tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                    "sample_E", "sample_F", "sample_G", "sample_H",
#'                    "sample_I", "sample_J", "sample_K", "sample_L",
#'                    "sample_M", "sample_N", "sample_O", "sample_P"),
#'                    identifier = "AS_events")
#' ts_sub(res2, Ib = 2, subtissues=c("sample_A", "sample_B"))
#' ts_sub(res2, Ib = 3, specificity = c(0.8, 1),
#'        subtissues = c("sample_A", "sample_B", "sample_C"))
#' ts_sub(res2, Ib = 15)
ts_sub <- function(lst,
                   Ib = 1,
                   specificity = c(0, 1),
                   ts_method = "Tau",
                   subtissues = NULL,
                   sample_mean = c(-Inf, Inf),
                   sample_mean_sd = Inf,
                   sample_max = c(-Inf, Inf),
                   sample_range = 0,
                   rep_diff_max = Inf,
                   rep_diff_sd = Inf,
                   rep_diff_mean = Inf) {
  if (!is.list(lst) | !(length(lst) == 5)) {
    stop("lst maybe fault input data!")
  }

  if (!is.numeric(Ib)) {
    stop("Ib should be numeric!")
  } else {
    Ib <- as.integer(Ib)
  }

  if (!is.numeric(sample_mean) | !is.numeric(sample_max) | !is.numeric(sample_range) | !is.numeric(rep_diff_max) | !is.numeric(sample_mean_sd) | !is.numeric(rep_diff_sd) | !is.numeric(rep_diff_mean)) {
    stop("sample_mean, sample_max, sample_range, sample_mean_sd, rep_diff_max, rep_diff_sd and rep_diff_mean should be numeric!")
  }

  if ((length(sample_mean) != 2) | (length(sample_max) != 2) | (length(sample_range) != 1) | (length(rep_diff_max) != 1) | (length(sample_mean_sd) != 1) | (length(rep_diff_sd) != 1) | (length(rep_diff_mean) != 1)) {
    stop("sample_mean, sample_max, sample_range, sample_mean_sd, rep_diff_max, rep_diff_sd and rep_diff_mean should be specified length!")
  }

  meths <- c("Tau", "Gini", "Tsi", "Counts", "Ee", "Pem", "Hg", "Z", "Spm")
  if ((length(ts_method) != 1) | !(ts_method %in% meths)) {
    stop("spe_meth error!")
  }

  df1 <- lst$raw
  df2 <- lst$rank
  df3 <- lst$bin
  df4 <- lst$diff
  df5 <- lst$origin
  num <- which(colnames(df1) == "Tau") - 1
  tissues <- colnames(df1)[1:num]
  #diff_tissue <- setdiff(tissues, subtissues)

  ## Ib
  if (is.null(subtissues)) {
    bin_genes <- rownames(df3[df3$Ib == Ib,])
  } else if (length(subtissues) != Ib) {
    stop("subtissues should has the length of Ib or assign NULL!")
  } else {
    bin_genes <- rownames(df3[df3$Ib == Ib,])
    bin_genes <- bin_genes[unname(apply(data.frame(df3[bin_genes, subtissues] == 1), 1, all))]
    #bin_gene1 <- unname(apply((df3[, subtissues] == rep(1, length(subtissues))), 1, all))
    #bin_gene2 <- unname(apply((df3[, diff_tissue] == rep(0, length(diff_tissue))), 1, all))
    #bin_genes <- rownames(df3[bin_gene1 & bin_gene2, ])
  }

  if (length(bin_genes) < 1) {
    message(paste("There is no genes of Ib =", Ib))
  }

  if (!all(subtissues %in% colnames(df3)[1:num])) {
    stop("subtissues error!")
  }

  if (!is.numeric(specificity) | (length(specificity) != 2)) {
    stop("specificity should be numeric!")
  } else {
    reg <- range(specificity)
    if (!(reg[1] >= 0) | (!reg[2] <= 1)) {
      stop("specificity error!")
    }
  }

  ## specificity
  spe_genes <- rownames(df1[which((df1[, ts_method] >= reg[1]) & (df1[, ts_method] <= reg[2])),])
  if (length(spe_genes) < 1) {
    message("There is no genes of specified specificity")
  }

  ## sample_mean = c(-Inf, Inf)
  reg <- range(sample_mean)
  sample_mean <- rownames(df1[which((df1[, "Mean"] >= reg[1]) & (df1[, "Mean"] <= reg[2])),])
  if (length(sample_mean) < 1) {
    message("There is no genes of specified sample_mean!")
  }

  ## sample_mean_sd = Inf
  sample_mean_sd <- rownames(df1[which(df1[, "sd"] <= sample_mean_sd),])
  if (length(sample_mean_sd) < 1) {
    message("There is no genes under specified sample_mean_sd!")
  }

  ## sample_max = c(-Inf, Inf)
  reg <- range(sample_max)
  sample_max <- rownames(df1[which((df1[, "Max"] >= reg[1]) & (df1[, "Max"] <= reg[2])),])
  if (length(sample_max) < 1) {
    message("There is no genes of specified sample_max!")
  }

  ## sample_range = -Inf
  sample_range <- rownames(df1[which(df1[, "Max_diff"] >= sample_range),])
  if (length(sample_range) < 1) {
    message("There is no genes of specified sample_range!")
  }

  ## rep_diff_max = Inf
  rep_diff_max <- rownames(df4[which(df4[, "Max_diff"] <= rep_diff_max),])
  if (length(rep_diff_max) < 1) {
    message("There is no genes of specified rep_diff_max!")
  }

  ## rep_diff_sd = Inf
  rep_diff_sd <- rownames(df4[which(df4[, "diff_sd"] <= rep_diff_sd),])
  if (length(rep_diff_sd) < 1) {
    message("There is no genes of specified rep_diff_sd!")
  }

  ## rep_diff_mean = Inf
  rep_diff_mean <- rownames(df4[which(df4[, "Mean_diff"] <= rep_diff_mean),])
  if (length(rep_diff_mean) < 1) {
    message("There is no genes of specified rep_diff_mean!")
  }

  int_genes <- Reduce(intersect, list(spe_genes, bin_genes, sample_mean, sample_max, sample_range, sample_mean_sd, rep_diff_max, rep_diff_sd, rep_diff_mean))

  if (length(int_genes) < 1) {
    new_lst <- list(raw = df1[int_genes, ],
                    rank = df2[int_genes, ],
                    bin = df3[int_genes, ],
                    diff = df4[int_genes,],
                    origin = df5[int_genes,]
    )
    message("There has no gene of the limitations!")
  } else {
    new_lst <- list(raw = df1[int_genes, ],
                    rank = df2[int_genes, ],
                    bin = df3[int_genes, ],
                    diff = df4[int_genes,],
                    origin = df5[int_genes,]
    )
    # if (is.null(subtissues)) {
    #   new_lst$bin$Type <- paste0(subtissues, collapse = ",")
    #   new_lst$rank$Type <- paste0(subtissues, collapse = ",")
    # } else {
      type <- apply(new_lst$bin[, 1:num], 1, function(x) {
        return(paste0(tissues[which(x == 1)], collapse = ","))
        })
      new_lst$bin$Type <- type
      new_lst$rank$Type <- type
    # }
    new.order <- order(-new_lst$raw$Tau, -new_lst$raw$Max_diff)
    new_lst <- lapply(new_lst, function(x) {return(x[new.order,])})
  }
  return(new_lst)
}


