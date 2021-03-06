% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{cmp_heatmap}
\alias{cmp_heatmap}
\title{Heatmap for compare rank between origin value (mean and replicates)}
\usage{
cmp_heatmap(lst, sort_fun = "de", trans = "none", main1 = "rank",
  main2 = "origin value", main3 = "origin value (Z-score)",
  cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE,
  ...)
}
\arguments{
\item{lst}{List of 5 data.frame. Generated from \code{ts_psi} or \code{ts_expr}.}

\item{sort_fun}{Heatmap row order. "in" (increase) or "de" (decrease).}

\item{trans}{Transformation of origion data for heatmap. "none" or "log2".}

\item{main1}{Title of rank plot.}

\item{main2}{Title of values (for calculate specificty) with replicates. For
expression, log2(CPM/FPKM/RPKM/TPM+1).}

\item{main3}{Title of origin values with Z-score.}

\item{cluster_cols}{refer to the same parameter in \code{pheatmap}.}

\item{cluster_rows}{refer to the same parameter in \code{pheatmap}.}

\item{show_rownames}{refer to the same parameter in \code{pheatmap}.}

\item{...}{parameters of \code{pheatmap} pkg.}
}
\description{
Wrapper function derived from \code{pheatmap}. Use the results of
\code{ts_sub} as Input.
}
\examples{
res <- ts_psi(demo_psi,
              tissues = c("sample_A", "sample_B", "sample_C",
                          "sample_D", "sample_E", "sample_F",
                          "sample_G", "sample_H", "sample_I",
                          "sample_J", "sample_K", "sample_L",
                          "sample_M", "sample_N", "sample_O",
                          "sample_P", "sample_Q"),
                          identifier = "AS_events")
sub_1 <- ts_sub(res, Ib = 1)
cmp_heatmap(sub_1,
            sort_fun = "de",
            main1 = "PSI rank",
            main2 = "PSI mean",
            main3 = "PSI (Z-score)")
}
