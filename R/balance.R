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
  # formulaBalance <- (group~variable+var1)
  # data <- dat
  # match_id = result$match_id

  #Check types of inputs (same function used for polymatch - amend with useless arguments)
  checkInputs(formulaMatch = formulaBalance, data = data, start = match_id,
              distance = "euclidean", iterate = T, niter_max = 50, verbose = T)

  #Check coherence of data (as above)
  resultCheckData <- checkData(formulaMatch = formulaBalance, data = data,
                               start = match_id)
  varGroup <- resultCheckData$varGroup
  varsBalance <- resultCheckData$varsMatch


  pairGroups <- combn(names(table(data[,varGroup])), 2)
  pairsGroupsText <- apply(pairGroups, FUN = function(x) {paste(x, collapse = "-")},2)

  #Generate a dataset to store the results of the balance: each variable has measure of balance
  # for each pair of groups
  dataBalance <- expand.grid(list(groups = pairsGroupsText,
                                  variable = varsBalance),
                             stringsAsFactors = F)

  dataBalance$type <- NA
  dataBalance$stdzDiffPre <- NA
  dataBalance$ratioVarsPre <- NA
  dataBalance$stdzDiffPost <- NA
  dataBalance$ratioVarsPost <- NA

  for(indexVar in 1:length(varsBalance)) {

    varBalance <- varsBalance[indexVar]

    if(class(data[,varBalance]) == "character") {
      data[,varBalance] <- factor(data[,varBalance])
    }

    typeVariableIter <- typeVariable(data[,varBalance])
    dataBalance$type[dataBalance$variable %in% varBalance] <- typeVariableIter

    for(indexPair in 1:ncol(pairGroups)) {

      selectionIter <- (dataBalance$groups %in% pairsGroupsText[indexPair] &
                          dataBalance$variable %in% varBalance)
      if(typeVariableIter == "continuous") {

        resultBalancePre <- balanceContVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                            varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiffPre[selectionIter] <- resultBalancePre$stdzDiff
        dataBalance$ratioVarsPre[selectionIter] <- resultBalancePre$ratioVars

        resultBalancePost <- balanceContVar(data = data, varBalance = varBalance, match_id = match_id,
                                        varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiffPost[selectionIter] <- resultBalancePost$stdzDiff
        dataBalance$ratioVarsPost[selectionIter] <- resultBalancePost$ratioVars

      }

      if(typeVariableIter == "binary") {

        resultBalancePre <- balanceBinVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                           varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiffPre[selectionIter] <- resultBalancePre$stdzDiff

        resultBalance <- balanceBinVar(data = data, varBalance = varBalance, match_id = match_id,
                                        varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiffPost[selectionIter] <- resultBalance$stdzDiff

      }

      if(typeVariableIter == "categorical") {

        resultBalancePre <- balanceCatVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                          varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiffPre[selectionIter] <- resultBalancePre$stdzDiff


        resultBalance <- balanceCatVar(data = data, varBalance = varBalance, match_id = match_id,
                                       varGroup = varGroup, pairGroups = pairGroups[,indexPair])
        dataBalance$stdzDiffPost[selectionIter] <- resultBalance$stdzDiff

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
plot.balanceCondOptMatch <- function(dataBalance) {

  #Data for standardized difference
  keepVarsStdzDiff <- c("groups","variable","stdzDiffPre","stdzDiffPost")
  dataBalanceStdzDiff <- tidyr::gather(dataBalance[,keepVarsStdzDiff],
                                       key = pre_post,
                                       value = stdzDiff, - variable, - groups)
  dataBalanceStdzDiff$pre_post <- factor(dataBalanceStdzDiff$pre_post,
                                         levels = c("stdzDiffPost","stdzDiffPre"),
                                         labels = c("Post","Pre"))
  dataBalanceStdzDiff$variable <- factor(dataBalanceStdzDiff$variable,
                                         levels = unique(dataBalance$variable))
  #Data for ratio of variances
  keepVarsRatioVars <- c("groups","variable","ratioVarsPre","ratioVarsPost")
  dataBalanceRatioVars <- tidyr::gather(dataBalance[dataBalance$type == "continuous",keepVarsRatioVars],
                                       key = pre_post,
                                       value = ratioVars, - variable, - groups)
  dataBalanceRatioVars$pre_post <- factor(dataBalanceRatioVars$pre_post,
                                         levels = c("ratioVarsPost","ratioVarsPre"),
                                         labels = c("Post","Pre"))
  dataBalanceRatioVars$variable <- factor(dataBalanceRatioVars$variable,
                                          levels = unique(resultBalance$variable[dataBalance$type == "continuous"]))



  plotStdzDiff <- ggplot2::ggplot(data = dataBalanceStdzDiff) +
    ggplot2::geom_boxplot(ggplot2::aes(x = pre_post,
                                       y = stdzDiff,
                                       colour= pre_post)) +
    ggplot2::geom_jitter(ggplot2::aes(x = pre_post,
                                       y = stdzDiff,
                                       colour= pre_post), size = 2) +
    ggplot2::facet_wrap(~variable, dir = "v", strip.position = "left", ncol = 1)  +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title="Standardized Differences Before and After Matching") +
    ggplot2::ylab("Standardized Differences Before and After Matching") +
    ggplot2::geom_hline(yintercept=0,  colour = "red", size = 1.2, linetype = "dashed") +
    ggplot2::guides(fill=FALSE, colour = FALSE)


  plotRatioVars <- ggplot2::ggplot(data = dataBalanceRatioVars) +
    ggplot2::geom_boxplot(ggplot2::aes(x = pre_post,
                                       y = ratioVars,
                                       colour= pre_post))  +
    ggplot2::geom_jitter(ggplot2::aes(x = pre_post,
                                      y = ratioVars,
                                      colour= pre_post), size = 2) +
    ggplot2::facet_wrap(~variable, dir = "v", strip.position = "left", ncol = 1)  +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title="Ratio of Variances Before and After Matching") +
    ggplot2::ylab("Ratio of Variances") +
    ggplot2::geom_hline(yintercept=1,  colour = "red", size = 1.2, linetype = "dashed") +
    ggplot2::guides(fill=FALSE, colour = FALSE)


  #Ideas: jittered point, jittered points on the same row,
  # non-jittered points, jittered points on the same row

  gridExtra::grid.arrange(plotStdzDiff, plotRatioVars, ncol = 2)

  return(list(plotStdzDiff=plotStdzDiff,
              plotRatioVars = plotRatioVars))
}
