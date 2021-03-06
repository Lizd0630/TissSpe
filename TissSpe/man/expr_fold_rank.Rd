% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary_pattern.R
\name{expr_fold_rank}
\alias{expr_fold_rank}
\title{Ranking gene expression according to fold-change-intervals}
\usage{
expr_fold_rank(df, n = 12, min = 0.5, medium = log2(3/2), step = 1)
}
\arguments{
\item{df}{data.frame.}

\item{n}{number * 2 of medium intervals to generated (all \code{2n}+2).
Default 12. Values greater than n set to highest rank.}

\item{min}{expression under \code{min} will be set to 0. Default 0.5.}

\item{medium}{the medium interval of neighbour 2 steps. Default log2(3/2).}

\item{step}{interval width. Default 1, log2 fold change: log2(2) = 1.}
}
\value{
data.frame contain ranks from \code{df}.
}
\description{
Function to cut data.frame range into \code{2n+2} fold-change-interval
points (1,2,3,4,6,8,12,16,24,32,...) and , maximal point and minimal
point, then rank them for their value. This function is design for
processed data, like log2 transformed gene expression
(0,1,1+log(3/2),2,2+log2(3/2)...)value(RPKM/FPKM/TPM).
}
