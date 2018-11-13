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
#' For example, if there are four groups with labels "A","B","C" and "D", \code{start="D-B-A-C"} generates the starting sample
#' by matching groups "D" and "B", then units from "A" to the "D"-"B"pairs, then units from "C" to the "D"-"B"-"A" triplets.
#' Third, users can provide a matched set as starting condition and the algorithm will explore possible reductions in the total
#' distance. In this case, \code{start} must be a vector with length equal to the number of rows of \code{data} and
#' matched subjects must be flagged with the same value. See section 'Details' for further information.
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
#'
#' The steps of the algorithm are described with the following example. Consider a 4-group design with
#' groups labels "A", "B", "C" and "D". The algorithm requires a set of quadruplets as starting point. The argument \code{start} defines the approach to be used to
#' generate such a starting point. \code{polymatch} generates the starting point by sequentially using optimal two-group matching.
#' In the default setting (\code{start="small.to.large"}), the steps are: 1) optimally match the two smallest
#' groups; 2) optimally match the third smallest group to the pairs generated in the first step; 3) optimally match the last group
#' to the triplets generated in the second step. Notably, we can use the optimal two-group algorithm in steps 2) and 3) because they are
#' two-dimensional problems: the elements of one group on one hand, fixed matched sets on the other hand. The order of the
#' groups to be considered when generating the estarting point can be user-specified (e.g., \code{start="D-B-A-C"}).
#' In alternative, the user can provide a matched set that will be used as starting point.
#'
#' Given the starting matched set, the algorithm iteratively explores possible reductions in the total distance (if \code{iterate="TRUE"}).
#' The algorithm sequentially relaxes the connection to each group and rematches units for that group:
#' i) rematch "B-C-D" triplets within the quadruplets to units in group "A";
#' ii) rematch "A-C-D" triplets within the quadruplets to units in group "B";
#' iii) rematch "A-B-D" triplets within the quadruplets to units in group "C";
#' iv) rematch "A-B-C" triplets within the quadruplets to units in group "D".
#' If none of the sets of quadruplets generated in i)-iv) has smaller total distance than the starting point, the algorihm stops.
#' Otherwise, the set of quadruplets with smallest distance is seleceted and the process iterated, until no reduction in the total
#' distance is found or the number of maximum iterations is reached (\code{niter_max=50} by default).
#'
#' The total distance is defined as the sum of all the within-matched-set distances. The within-matched-set distance is defined as the
#' sum of the pairwise distances between pairs of units in the matched set. The type of distance is specified with the \code{distance}
#' argument. The current implementation supports Euclidean (\code{distance="euclidean"}) and Mahalanobis (\code{distance="mahalanobis"})
#' distances. In particular, for the Mahalanobis distance, the covariance matrix is defined only once on the full dataset.
#'
#' @examples
#' plot(1, 1)
#'
#' @export
polymatch <- function(formulaMatch, start = "small.to.large", data, distance = "euclidean", iterate = T, niter_max = 50, verbose = T) {

  #Debug/devel:
  #------------
  # source("C:/Users/natt03/Documents/R/temp.R")
  # formulaMatch <- (groupF~variable)
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
  vectorSchemeStart <- resultCheckData$vectorSchemeStart

  #Select variables of interest
  data <- data[,c(varGroup,varsMatch)]

  #Make grouping variable as character
  data[,varGroup] <- as.character(data[,varGroup])

  #Number of groups
  numGroups <- length(table(data[,varGroup]))

  #Define an id to sort observations in order provided
  data$idUnits <- 1:nrow(data)

  if(distance == "mahalanobis") {
    Sigma <- stats::cov(as.matrix(data[,varsMatch]))
  } else {
    Sigma <- NULL
  }

  # Say some stuff
  if(verbose==T) {

    cat("Conditional optimal matching algorithm\n")
    cat("Number of observations: ",nrow(data),"\n")
    cat("Number of groups: ",numGroups,"\n")
    #cat("\n")
  }

  #First: generate and evaluate starting point
  #-------------------------------------------

  #IF: starting matched dataset provided
  if( is.null(vectorSchemeStart)) {

    data$match_id <- start

    #Evaluate total distance at starting condition
    resultEvaluation <- evaluateMatching(data, "match_id", varsMatch,
                                         distance, Sigma)
    total_distance_start <- resultEvaluation$total_distance


    tabGroup <- table(data[,varGroup])
    vectorSchemeIter <- names(sort(tabGroup))

  #ELSE: starting matched dataset must be constructed
  } else {

    #First step: select units in first group
    data1 <- data[data[,varGroup] %in% vectorSchemeStart[1], ]

    #Each unit is a matched set with 1 element
    data1$indexMatch1 <- NA
    data1$indexMatch1[data1[,varGroup] %in% vectorSchemeStart[1]] <- 1:nrow(data1)

    for(i in 2:numGroups) {

      #Next step: select units from the next group
      data2 <- data[data[,varGroup] %in% vectorSchemeStart[i], ]

      #Each unit is a matched set with 1 element
      data2$indexMatch2 <- NA
      data2$indexMatch2[data2[,varGroup] %in% vectorSchemeStart[i]] <- 1:nrow(data2)

      #Before combining the dataset, variables need to be the same
      data2$indexMatch1 <- NA
      data1$indexMatch2 <- NA

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
    data$indexMatch1 <- NULL

    total_distance_start <- resultStep$total_distance
    vectorSchemeIter <- vectorSchemeStart

  }

  # Say some stuff
  if(verbose==T) {
    cat("Total distance of starting matched sample: ", sprintf("%.3f",total_distance_start),"\n")
    #cat("\n")
  }

  #Second: iterations
  #-------------------

  best_total_distance <- total_distance_start
  best_match_id <- data$match_id

  if(iterate) {

    for(niter in 1:niter_max) {

      data$match_id <- best_match_id

      improvementInIteration <- F

      for(iGroupStepIter2 in 1:length(vectorSchemeIter)) {

        groupStepIter2 <- vectorSchemeIter[iGroupStepIter2]
        groupsStepIter1 <- setdiff(vectorSchemeIter, groupStepIter2)

        #Relax connection to groupStepIter2
        data$indexMatchIter1 <- data$match_id
        data$indexMatchIter1[data[,varGroup] %in% groupStepIter2] <- NA

        data$indexMatchIter2 <- NA
        data$indexMatchIter2[data[,varGroup] %in% groupStepIter2] <- 1:sum(data[,varGroup] %in% groupStepIter2)

        #Rematch groupsStepIter1 to groupStepIter2
        resultIter <- condOptMatching(data = data,
                                      varIndexMatch1 = "indexMatchIter1",
                                      varIndexMatch2 = "indexMatchIter2",
                                      varsMatch = varsMatch,
                                      varGroup = varGroup,
                                      distance = distance,
                                      Sigma = Sigma)

        if(resultIter$total_distance < best_total_distance) {

          best_match_id <- resultIter$match_id
          best_total_distance <- resultIter$total_distance
          improvementInIteration <- T
        }

      }

      if(improvementInIteration == F) {
        break;
      }

      # Say some stuff
      if(verbose==T) {
        cat("Ended iteration ", niter, " - total distance:", sprintf("%.3f",best_total_distance),"\n")
        #cat("\n")
      }


    }

  } else {

    total_distance <- total_distance_start
    niter <- 1

  }

  if(verbose == T) {
    cat("End \n")
    cat("Number of iterations: ", niter, ", total distance:", sprintf("%.3f",best_total_distance),"\n" )
  }

  if(iterate == T & niter>=niter_max) {
    warning("The algorithm reached the maximum number of iterations--you can increase it with the argument 'niter_max'")
  }

  dataToOutput <- data[order(data$idUnits), c(varGroup,varsMatch)]

  return(list(match_id = best_match_id[order(data$idUnits)],
              total_distance = best_total_distance,
              niter = niter,
              total_distance_start = total_distance_start))
}
