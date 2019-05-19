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
#' @param subtissues Charactor, vector. Of binary pattern, only tissues in
#' \code{subtissues} have value 1, that is specific-overexpression.
#' @return a list with 3 data.frames. Subset of \code{lst}.
#' @export
#' @examples
#' res1 <- ts_expr(demo_tpm, n = 11,
#'         tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                    "sample_E", "sample_F", "sample_G", "sample_H",
#'                    "sample_I", "sample_J", "sample_K", "sample_L",
#'                    "sample_M", "sample_N", "sample_O", "sample_P"),
#'                    identifier = "gene_id")
#' ts_sub(res1, Ib=2, subtissues=c("sample_A", "sample_B"))
#' ts_sub(res1, Ib=3, specificity = c(0, 1),
#'        subtissues = c("sample_A", "sample_B", "sample_C"))
#' res2 <- ts_psi(demo_psi,
#'         tissues = c("sample_A", "sample_B", "sample_C", "sample_D",
#'                    "sample_E", "sample_F", "sample_G", "sample_H",
#'                    "sample_I", "sample_J", "sample_K", "sample_L",
#'                    "sample_M", "sample_N", "sample_O", "sample_P"),
#'                    identifier = "gene_id")
#' ts_sub(res2, Ib=2, subtissues=c("sample_A", "sample_B"))
#' ts_sub(res2, Ib=3, specificity = c(0.8, 1),
#'        subtissues = c("sample_A", "sample_B", "sample_C"))
ts_sub <- function(lst,
                   Ib = 1,
                   specificity = c(0, 1),
                   ts_method = "Tau",
                   subtissues) {
  if (!is.list(lst) | !(length(names(lst)) == 3)) {
    stop("lst maybe fault input data!")
  }

  reg <- range(specificity)
  if (!(is.numeric(specificity) & (length(reg) == 2) & (reg[1] >= 0) & (reg[2] <= 1))) {
    stop("specificity error!")
  }

  meths <- c("Tau", "Gini", "Tsi", "Counts", "Ee", "Pem", "Hg", "Z", "Spm")
  if ((length(ts_method) != 1) | !(ts_method %in% meths)) {
    stop("spe_meth error!")
  }

  if (length(subtissues) != Ib) {
    stop("subtissues should has the length of Ib!")
  }

  df1 <- lst$raw
  df2 <- lst$rank
  df3 <- lst$bin
  num <- ncol(df1) - 11
  tissues <- colnames(df1)[1:num]
  #diff_tissue <- setdiff(tissues, subtissues)

  if (!all(subtissues %in% colnames(df3)[1:num])) {
    stop("subtissues error!")
  }

  reg <- range(specificity)
  spe_genes <- rownames(df1[which((df1[, ts_method] >= reg[1]) & (df1[, ts_method] <= reg[2])),])
  if (length(spe_genes) < 1) {
    stop(paste("There is no genes of specificity", specificity))
  }

  bin_genes <- rownames(df3[df3$Ib == Ib,])
  bin_genes <- bin_genes[unname(apply((df3[bin_genes, subtissues] == rep(1, length(subtissues))), 1, all))]
  #bin_gene1 <- unname(apply((df3[, subtissues] == rep(1, length(subtissues))), 1, all))
  #bin_gene2 <- unname(apply((df3[, diff_tissue] == rep(0, length(diff_tissue))), 1, all))
  #bin_genes <- rownames(df3[bin_gene1 & bin_gene2, ])
  if (length(bin_genes) < 1) {
    stop(paste("There is no genes of Ib =", Ib))
  }

  int_genes <- intersect(spe_genes, bin_genes)

  new_lst <- list(raw = df1[int_genes, ],
                  rank = df2[int_genes, ],
                  bin = df3[int_genes, ])

  return(new_lst)
}


