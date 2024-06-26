#' Evaluating the Balance of Covariates After Matching
#'
#' The function \code{balance} computes the standardized mean differences and the ratio of the variances among treatment groups,
#' before and after matching. The function computes the two measures of balance for each pair of treatment groups.
#'
#' @param formulaBalance Formula with form \code{group ~ x_1 + ... + x_p}. \code{group} is the variable
#' identifying the treatment groups/exposures. The balance is evaluated for the covariates \code{x_1},...,\code{x_p}.
#' Numeric and integer variables are treated as continuous. Factor variables are treated as categorical.
#' Factor variables with two levels are treated as binary.
#' @param match_id Vector identifying the matched sets---matched units must have the same identifier. It is generated by
#' \code{\link{polymatch}}.
#' @param data The \code{data.frame} object with the data.
#' @param weights_before Optional vector of weights of the observations to be considered in the unmatched dataset. To compute the
#' unweighted standardized mean differences, set \code{weights_before} to NULL (default).
#' @param weights_after Vector of weights for the matched dataset. Set it to NULL (default) to compute the
#' unweighted standardized mean differences.
#'
#' @return A \code{data.frame} containing the standardized differences and ratios of the variances (only for continuous
#' variables) for each pair of treatment groups. A graphical representation of the results can be generated with
#' \code{\link{plotBalance}}.
#'
#' @seealso \code{\link{polymatch}} to generate matched samples and \code{\link{plotBalance}} to
#' graphically represent the indicators of balance.
#'
#' @examples
#' #Generate a datasets with group indicator and four variables:
#' #- var1, continuous, sampled from normal distributions;
#' #- var2, continuous, sampled from beta distributions;
#' #- var3, categorical with 4 levels;
#' #- var4, binary.
#' set.seed(1234567)
#' dat <- data.frame(group = c(rep("A",20),rep("B",60),rep("C",60)),
#'                   var1 = c(rnorm(20,mean=0,sd=1),
#'                            rnorm(60,mean=1,sd=2),
#'                            rnorm(60,mean=-1,sd=2)),
#'                   var2 = c(rbeta(20,shape1=1,shape2=1),
#'                            rbeta(60,shape1=2,shape2=1),
#'                            rbeta(60,shape1=1,shape2=2)),
#'                   var3 = factor(c(rbinom(20,size=3,prob=.4),
#'                                   rbinom(60,size=3,prob=.5),
#'                                   rbinom(60,size=3,prob=.3))),
#'                   var4 = factor(c(rbinom(20,size=1,prob=.5),
#'                                   rbinom(60,size=1,prob=.3),
#'                                   rbinom(60,size=1,prob=.7))))
#'
#' #Match on propensity score
#' #-------------------------
#'
#' #With multiple groups, need a multinomial model for the PS
#' library(VGAM)
#' psModel <- vglm(group ~ var1 + var2 + var3 + var4,
#'                 family=multinomial, data=dat)
#' #Estimated probabilities - 3 for each unit: P(group=A), P(group=B), P(group=C)
#' probsPS <- predict(psModel, type = "response")
#' dat$probA <- probsPS[,"A"]
#' dat$probB <- probsPS[,"B"]
#' dat$probC <- probsPS[,"C"]
#' #Estimated logits - 2 for each unit: log(P(group=A)/P(group=C)), log(P(group=B)/P(group=C))
#' logitPS <- predict(psModel, type = "link")
#' dat$logit_AvsC <- logitPS[,1]
#' dat$logit_BvsC <- logitPS[,2]
#'
#' #Match on logits of PS
#' resultPs <- polymatch(group ~ logit_AvsC + logit_BvsC, data = dat,
#'                     distance = "euclidean")
#' dat$match_id_ps <- resultPs$match_id
#'
#' #Evaluate balance in covariates
#' tabBalancePs <- balance(group ~ var1 + var2 + var3 + var4,
#'                         match_id = dat$match_id_ps, data = dat)
#' tabBalancePs
#'
#' #You can also represent the standardized mean differences with 'plotBalance'
#' #plotBalance(tabBalancePs, ratioVariances = TRUE)
#'
#' @export
balance <- function(formulaBalance, match_id, data, weights_before = NULL, weights_after = NULL) {

  #browser()
  #Debug/devel:
  #------------
  # formulaBalance <- (group~var1)
  # data <- dat
  # match_id = result$match_id

  #Check types of inputs (same function used for polymatch - amend with useless arguments)
  checkInputs(formulaMatch = formulaBalance, start = match_id, data = data,
              distance = "euclidean", exactMatch = NULL, iterate = TRUE, niter_max = 50, verbose = TRUE, vectorK = NULL)

  #Check coherence of data (as above)
  resultCheckData <- checkData(formulaMatch = formulaBalance, start = match_id,
                               data = data, exactMatch = NULL, vectorK = NULL, checkOnePerGroup = FALSE)
  varGroup <- resultCheckData$varGroup
  varsBalance <- resultCheckData$varsMatch
  vectorK <- resultCheckData$vectorK

  #Add weights to the dataset
  #Before matching
  if (!is.null(weights_before)) {
    data$weights_before <- weights_before
  } else {
    data$weights_before <- rep(1, nrow(data))
  }
  varWeightsBefore <- "weights_before"
  #After matching
  if (!is.null(weights_after)) {
    data$weights_after <- weights_after
  } else {
    data$weights_after <- as.vector(1/vectorK)[match(data[,varGroup], names(vectorK))]
  }
  varWeightsAfter <- "weights_after"

  #Define all pairs of groups
  pairGroups <- utils::combn(names(table(data[,varGroup])), 2)
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
                                           varGroup = varGroup, pairGroups = pairGroups[,indexPair],
                                           varWeights = varWeightsBefore)
        dataBalance$stdzDiffPre[selectionIter] <- resultBalancePre$stdzDiff
        dataBalance$ratioVarsPre[selectionIter] <- resultBalancePre$ratioVars

        resultBalancePost <- balanceContVar(data = data, varBalance = varBalance, match_id = match_id,
                                            varGroup = varGroup, pairGroups = pairGroups[,indexPair],
                                            varWeights = varWeightsAfter)
        dataBalance$stdzDiffPost[selectionIter] <- resultBalancePost$stdzDiff
        dataBalance$ratioVarsPost[selectionIter] <- resultBalancePost$ratioVars

      }

      if(typeVariableIter == "binary") {

        resultBalancePre <- balanceBinVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                          varGroup = varGroup, pairGroups = pairGroups[,indexPair],
                                          varWeights = varWeightsBefore)
        dataBalance$stdzDiffPre[selectionIter] <- resultBalancePre$stdzDiff

        resultBalance <- balanceBinVar(data = data, varBalance = varBalance, match_id = match_id,
                                       varGroup = varGroup, pairGroups = pairGroups[,indexPair],
                                       varWeights = varWeightsAfter)
        dataBalance$stdzDiffPost[selectionIter] <- resultBalance$stdzDiff

      }

      if(typeVariableIter == "categorical") {

        resultBalancePre <- balanceCatVar(data = data, varBalance = varBalance, match_id = rep(1,nrow(data)),
                                          varGroup = varGroup, pairGroups = pairGroups[,indexPair],
                                          varWeights = varWeightsBefore)
        dataBalance$stdzDiffPre[selectionIter] <- resultBalancePre$stdzDiff


        resultBalance <- balanceCatVar(data = data, varBalance = varBalance, match_id = match_id,
                                       varGroup = varGroup, pairGroups = pairGroups[,indexPair],
                                       varWeights = varWeightsAfter)
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
#' @param dataBalance the output of \code{\link{balance}}.
#' @param ratioVariances Boolean. If \code{TRUE}, the generated plot contains two panels:
#' one for the standardized differences and one for the ratios of the variances. If \code{FALSE}
#' (the default), only the standardized differences are represented.
#' @param boxplots Boolean. If \code{TRUE} (default), boxplots are added to the plot, to show the
#' distribution of the standardized differences and ratios of the variances.
#'
#' @return If at least one of the covariates is continuous and \code{ratioVariances=TRUE},
#' the function generates a plot with two panels: one for the
#' standardized differences and one for the ratio of the variances (only for the continous variables).
#' If either all the covariates are categorical/binary or \code{ratioVariances=FALSE} (or both),
#' only the plot with the standardized differences is generated.
#' The function also returns a list with the \code{ggplot2} objects corresponding to the generated plot(s).
#'
#' @seealso \code{\link{polymatch}} to generate matched samples and \code{\link{balance}} to compute
#' the indicators of balance.

#' @examples
#' #See examples of function 'balance'
#'
#' @export
plotBalance <- function(dataBalance, ratioVariances = FALSE, boxplots = TRUE) {

  #Data for standardized difference
  keepVarsStdzDiff <- c("groups","variable","stdzDiffPre","stdzDiffPost")
  dataBalanceStdzDiff <- tidyr::gather(dataBalance[,keepVarsStdzDiff],
                                       key = "pre_post",
                                       value = "stdzDiff", - "variable", - "groups")
  dataBalanceStdzDiff$pre_post <- factor(dataBalanceStdzDiff$pre_post,
                                         levels = c("stdzDiffPost","stdzDiffPre"),
                                         labels = c("Post","Pre"))
  dataBalanceStdzDiff$variable <- factor(dataBalanceStdzDiff$variable,
                                         levels = unique(dataBalance$variable))

  #Data for ratio of variances
  if(ratioVariances == TRUE) {

    keepVarsRatioVars <- c("groups","variable","ratioVarsPre","ratioVarsPost")
    dataBalanceRatioVars <- tidyr::gather(dataBalance[dataBalance$type == "continuous",keepVarsRatioVars],
                                         key = "pre_post",
                                         value = "ratioVars", - "variable", - "groups")
    dataBalanceRatioVars$pre_post <- factor(dataBalanceRatioVars$pre_post,
                                           levels = c("ratioVarsPost","ratioVarsPre"),
                                           labels = c("Post","Pre"))
    dataBalanceRatioVars$variable <- factor(dataBalanceRatioVars$variable,
                                            levels = unique(dataBalance$variable[dataBalance$type == "continuous"]))
  } else {
    #Empty data frame with 0 variables if ratioVariances = F
    dataBalanceRatioVars <- data.frame(var=numeric(0))
  }


  #How many columns? If also the ratio of the variances: 1 column
  if(nrow(dataBalanceRatioVars)>0) {

    numberColumnsStdzDiff <- 1


    } else {

      #If only standardized differnces: 1 column if not many variables.
      if(length(unique(dataBalanceStdzDiff$variable))<=9) {

        numberColumnsStdzDiff <- 1

      } else {
        #If only standardized differnces but many variables: 2 columns.
        numberColumnsStdzDiff <-   2

      }

  }

  plotStdzDiff <- ggplot2::ggplot(data = dataBalanceStdzDiff) +
    ggplot2::aes_string(x = "pre_post",
                        y = "stdzDiff",
                        colour= "pre_post") +
    ggplot2::geom_jitter(size = 2) +
    ggplot2::facet_wrap(~variable, dir = "v", strip.position = "left", ncol = numberColumnsStdzDiff)  +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::labs(title="Standardized Differences Before and After Matching") +
    ggplot2::ylab("Standardized Differences") +
    ggplot2::geom_hline(yintercept=0,  colour = "red", size = 1.2, linetype = "dashed") +
    ggplot2::guides(fill=FALSE, colour = FALSE)

  if(boxplots==T) {
    plotStdzDiff <- plotStdzDiff +
      ggplot2::geom_boxplot(ggplot2::aes_string(x = "pre_post",
                                                  y = "stdzDiff",
                                                  colour= "pre_post"),
                            outlier.shape = NA)
    }

  #Plot of ratio of variances only if there is at least one continuous variable
  if(nrow(dataBalanceRatioVars)>0) {

    plotRatioVars <- ggplot2::ggplot(data = dataBalanceRatioVars) +
      ggplot2::geom_jitter(ggplot2::aes_string(x = "pre_post",
                                        y = "ratioVars",
                                        colour= "pre_post"), size = 2) +
      ggplot2::facet_wrap(~variable, dir = "v", strip.position = "left", ncol = 1)  +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::labs(title="Ratio of Variances Before and After Matching") +
      ggplot2::ylab("Ratio of Variances") +
      ggplot2::geom_hline(yintercept=1,  colour = "red", size = 1.2, linetype = "dashed") +
      ggplot2::guides(fill=FALSE, colour = FALSE)

    if(boxplots==T) {
      plotRatioVars <- plotRatioVars +
        ggplot2::geom_boxplot(ggplot2::aes_string(x = "pre_post",
                                                  y = "ratioVars",
                                                  colour= "pre_post"),
                              outlier.shape = NA)
    }

  } else {

    plotRatioVars <- NULL

  }

  if(nrow(dataBalanceRatioVars)>0) {

    gridExtra::grid.arrange(plotStdzDiff, plotRatioVars, ncol = 2)

  } else {

    print(plotStdzDiff)

  }


  return(list(plotStdzDiff=plotStdzDiff,
              plotRatioVars = plotRatioVars))
}
