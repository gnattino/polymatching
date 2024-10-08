% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matching_functions.R
\name{condOptMatching}
\alias{condOptMatching}
\title{Function implementing conditional optimal matching.}
\usage{
condOptMatching(
  data,
  varIndexMatch1,
  varIndexMatch2,
  varsMatch,
  varGroup,
  distance,
  Sigma,
  varsExactMatch,
  k,
  dat_stdzDistances
)
}
\description{
Optimally matches two sets of IDs of matched units--varIndexMatch1 and varIndexMatch2--with respect
to matching variables varsMatch.
}
\keyword{internal}
