#' A psi demo of AS events in 17 tissues
#'
#' A dataset containing the alternative splicing events and their
#' percentage-of-splicing-in.
#'
#' @format A data.frame with 5000 rows and 35variables:
#' \describe{
#'   \item{AS_events}{Alternative splicing has same start}
#'   \item{sample_A_1}{sample_A, replicate 1}
#'   \item{sample_A_2}{sample_A, replicate 2}
#'   ...
#' }
"demo_psi"

#' A gene expression demo of TPM in 16 tissues
#'
#' A dataset containing the gene names and their expression in tissues of TPM.
#'
#' @format A data.frame with 5000 rows and 33 variables:
#' \describe{
#'   \item{gene_id}{ENSEMBL gene id}
#'   \item{sample_A_1}{sample_A, replicate 1}
#'   \item{sample_A_2}{sample_A, replicate 2}
#'   ...
#' }
"demo_tpm"
