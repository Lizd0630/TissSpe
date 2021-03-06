% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spe_meth.R
\name{ts_Z}
\alias{ts_Z}
\title{Z-score}
\usage{
ts_Z(x)
}
\arguments{
\item{x}{data.frame.}
}
\description{
Function require a vector with expression of one gene in different tissues.
If expression for one tissue is not known, gene specificity for this gene is
NA. Transformation: x/n-1/sqrt(N).
}
\examples{
\dontrun{
ts_Z(c(0,1,2,3,2,8))
}
}
