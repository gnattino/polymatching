#' Polymatching
#'
#' \code{polymatch} is the main function to generate matched samples in designs with 3, 4, 5 and 6 groups.
#'
#' @param formulaMatch Formula with form \code{group ~ x_1 + ... + x_p}, where \code{group} is the name of the variable
#' identifying the treatment groups/exposures and \code{x_1},...,\code{x_p} are the matching variables.
#' @param data The \code{data.frame} object with the data. It must contain all the variables specified in \code{formulaMatch}.
#' @param distance String specifying whether the distance between pairs of observations should be computed with the
#' Euclidean (\code{"euclidean"}, default) or Mahalanobis (\code{"mahalanobis"}) distance. See section 'Details' for further information.
#' @param start An object specifying the starting point of the iterative algorithm. Three types of inputs are accepted.
#' First, \code{start="small.to.large"} (default) generates the first set of matched sets by matching groups
#' from the smallest to the largest. Second, users can specify the order to be used to match groups for the starting sample.
#' For example, if there are four groups with labels "a","b","c" and "d", \code{start="d-b-a-c"} generates the starting sample
#' by matching groups "d" and "b", than units from "a" to the "d"-"b"pairs, then units from "c" to the "d"-"b"-"a" triplets.
#' Third, users can provide a matched set as starting condition and the algorithm will explore possible reductions in the total
#' distance. In this case, \code{start} must be a vector with length equal to the number of rows of \code{data} and
#' matched subjects must be flagged with the same value.
#' @param iterate Boolean specifying whether iterations should be done (\code{iterate=TRUE}, default) or not (\code{iterate=FALSE}).
#' @param niter_max Maximum number of iterations. Default is 50.
#' @param verbose Boolean: should text be printed in the console? Default is \code{TRUE}.

#' @return A list containing the following components:
#' \describe{
#'   \item{match_id}{A numeric vector identifying the matched sets--matched units have the same ID}
#'   \item{total_distance}{Total distance of the returned matched sample.}
#'   \item{total_distance_start}{Total distance at the starting point.}
#' }
#'
#' @details The function implements the conditionally optimal matching algorithm, which iteratively uses
#' two-group optimal matching to generate matched samples with small total distance. In the current implementation,
#' it is possible to generate matched samples with only one subject per group.
#' Describe distance within matched sets.
#'
#' @examples
#' plot(NA, NA)
#'
#' @export


polymatch <- function(formulaMatch, data, distance = "euclidean", start = "small.to.large", iterate = T, niter_max = 50, verbose = T) {

  #Debug/devel:
  #------------
  # source("C:/Users/natt03/Documents/R/temp.R")
  # formulaMatch <- (group~variable)
  # data <- generateData(c(30,10,40,20))
  # distance = "euclidean"
  # start = "small.to.large"
  # iterate = T
  # niter_max = 50
  # verbose = T

  #Check types of inputs
  checkInputs(formulaMatch, data, distance, start, iterate, niter_max, verbose)

  #Check coherence of data
  resultCheckData <- checkData(formulaMatch, data, start)
  varGroup <- resultCheckData$varGroup
  varsMatch <- resultCheckData$varsMatch
  vectorScheme <- resultCheckData$vectorScheme

  if( is.null(vectorScheme)) {
    #start from provided matched set
  } else {
    #generate  starting point
  }



}
