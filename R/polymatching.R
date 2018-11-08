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

  data <- data[,c(varGroup,varsMatch)]

  if(distance == "mahalanobis") {
    Sigma <- cov(dataStep[,varsMatch])
  } else {
    Sigma <- NULL
  }

  #IF: starting matched dataset provided
  if( is.null(vectorScheme)) {

    data$match_id <- start

    #Evaluate total distance at starting condition
    resultEvaluation <- evaluateMatching(data, "match_id", varsMatch,
                                         distance, Sigma)
    total_distance_start <- resultEvaluation$total_distance

  #ELSE: starting matched dataset must be constructed
  } else {

    numGroups <- length(vectorScheme)

    #First step: select units in first group
    data1 <- data[data[,varGroup] %in% vectorScheme[1], ]

    #Each unit is a matched set with 1 element
    data1$indexMatch1 <- NA
    data1$indexMatch1[data1[,varGroup] %in% vectorScheme[1]] <- 1:nrow(data1)

    for(i in 2:numGroups) {

      #Next step: select units from the next group
      data2 <- data[data[,varGroup] %in% vectorScheme[i], ]

      #Each unit is a matched set with 1 element
      data2$indexMatch2 <- NA
      data2$indexMatch2[data2[,varGroup] %in% vectorScheme[i]] <- 1:nrow(data2)

      #Append data of new group to previous data
      dataStep <- rbind(data1, data2)

      #Optimally match unit of data2 to matched sets in data1
      resultStep <- condOptMatching(data = dataStep,
                                    varIndexMatch1 = "indexMatch1",
                                    varIndexMatch2 = "indexMatch2",
                                    varsMatch = varsMatch,
                                    varGroup = varGroup,
                                    distance = distance,
                                    Sigma = Sigma)

      #Erase ids of previous matching step
      dataStep$indexMatch1 <- dataStep$indexMatch2 <- NULL

      #Id of new matching step.
      dataStep$indexMatch1 <- resultStep$match_id
      data1 <- dataStep

    }

    data <- dataStep
    data$match_id <- data$indexMatch1
    data$indexMatch1 < NULL

    total_distance_start <- resultStep$total_distance

  }



}
