% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covcpd_utils_function.R
\name{CountMap6}
\alias{CountMap6}
\title{Generate a scRNA-seq count matrix, allowing the input of mean expression levels and dropout probabilities.
See C0301, C0302 for function usages.}
\usage{
CountMap6(
  sigma,
  ngene,
  n,
  a1 = 3,
  b1 = 1,
  a20 = 2,
  b20 = 3,
  a30 = 1,
  b30 = 10,
  scmeanvec = NULL,
  pijvec = NULL,
  meanDiffPct = 0
)
}
\value{
a list of covariance matrices.
}
\description{
Generate a scRNA-seq count matrix, allowing the input of mean expression levels and dropout probabilities.
See C0301, C0302 for function usages.
}
