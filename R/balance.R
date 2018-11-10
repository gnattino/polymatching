#' Evaluation of Balance in Covariates After Matching
#'
#' \code{balance} .
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
balance <- function(formulaBalance, match_id, data) {

  #Debug/devel:
  #------------
  # formulaBalance <- (group~var1+var2)
  # data <- dat
  # match_id = resultM$match_id

  #Check types of inputs (same function used for polymatch - amend with useless arguments)
  checkInputs(formulaMatch = formulaBalance, data = data, start = match_id,
              distance = "euclidean", iterate = T, niter_max = 50, verbose = T)

  #Check coherence of data (as above)
  resultCheckData <- checkData(formulaMatch = formulaBalance, data = data,
                               start = match_id)
  varGroup <- resultCheckData$varGroup
  varsBalance <- resultCheckData$varsMatch

  dataBalance <- data.frame(variable = varsBalance,
                            type = NA,
                            stdzDiffPre = NA,
                            ratioVarsPre = NA,
                            stdzDiffPost = NA,
                            ratioVarsPost = NA,
                            stringsAsFactors = F)

  pairGroups <- combn(names(table(data[,varGroup])), 2)

  for(i in 1:lenght(varsBalance)) {

    varBalance <- dataBalance$variable[i]

    dataBalance$type[i] <- typeVariable(data[,varBalance])

    for(indexPair in 1:ncol(pairGroups)) {

      if(dataBalance$type[i] == "continuous") {

        resultBalancePre <- balanceContVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                            varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDifPre[i] <- resultBalancePre$stdzDiff
        dataBalance$ratioVarsPre[i] <- resultBalancePre$ratioVarsPost

        resultBalancePost <- balanceContVar(data = data, varBalance = varBalance, match_id = match_id,
                                        varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDifPost[i] <- resultBalancePost$stdzDiff
        dataBalance$ratioVarsPost[i] <- resultBalancePost$ratioVarsPost

      }

      if(dataBalance$type[i] == "binary") {

        resultBalancePre <- balanceBinVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                           varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDifPre[i] <- resultBalancePre$stdzDiff

        resultBalance <- balanceBinVar(data = data, varBalance = varBalance, match_id = match_id,
                                        varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiff[i] <- resultBalance$stdzDiff

      }

      if(dataBalance$type[i] == "categorical") {

        resultBalancePre <- balanceCatVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                          varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDifPre[i] <- resultBalancePre$stdzDiff


        resultBalance <- balanceCatVar(data = data, varBalance = varBalance, match_id = match_id,
                                       varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiff[i] <- resultBalance$stdzDiff

      }


    }

  }


  class(dataBalance) <- c(class(dataBalance), "balanceCondOptMatch")

  return(dataBalance)
}

#' Summary Plot of Balance in Covariates
#'
#' The function generates a plot summarizing the balance of the covariates.
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
plot.balanceCondOptMatch <- function() {

}
