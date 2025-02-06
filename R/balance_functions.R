#' Check type of variable
#' @keywords internal
typeVariable <- function(variable) {

  if(methods::is(variable, "numeric") | methods::is(variable, "integer")) {

    type <- "continuous"

  } else {

    if(methods::is(variable, "factor")) {

      if(length(levels(variable))<=2) {

        type <- "binary"

      } else {

        type <- "categorical"
      }

    } else {

      stop("Covariates must be numeric or factor variables")

    }

  }

  return(type)
}

#' Function to compute weighted variance
#' @keywords internal
weighted.variance <- function(x, w) {

  sum_w <- sum(w)
  sum_wsq <- sum(w^2)
  mean <- stats::weighted.mean(x = x, w = w)

  sum(w * (x - mean)^2) * (sum_w / (sum_w^2 - sum_wsq))

}


#' Balance indicators for continuous variables
#' @keywords internal
balanceContVar <- function(data, varBalance, match_id, varGroup, pairGroups, varWeights){

  varGroup_symbol <- as.symbol(varGroup)
  varBalance_symbol <- as.symbol(varBalance)
  varWeights_symbol <- as.symbol(varWeights)

  #For difference in means, matched data
  dataMatch <- data[data[,varGroup] %in% pairGroups & !is.na(match_id), ]
  infoMatch <- dataMatch %>%
    dplyr::group_by(!!varGroup_symbol) %>%
    dplyr::summarize(means = stats::weighted.mean(x = !!varBalance_symbol, w = !!varWeights_symbol),
                     vars = weighted.variance(x = !!varBalance_symbol, w = !!varWeights_symbol))

  #For variances, unmatched data
  dataUnm <- data[data[,varGroup] %in% pairGroups, ]
  infoUnm <- dataUnm %>%
    dplyr::group_by(!!varGroup_symbol) %>%
    dplyr::summarize(vars = weighted.variance(x = !!varBalance_symbol, w = !!varWeights_symbol))

  stdzDiff <- (infoMatch$means[infoMatch[,varGroup]==pairGroups[1]] - infoMatch$means[infoMatch[,varGroup]==pairGroups[2]])/(
                sqrt((infoUnm$vars[infoUnm[,varGroup]==pairGroups[1]] + infoUnm$vars[infoUnm[,varGroup]==pairGroups[2]])/2))

  #Return the ratio of the variances >= 1
  ratioVars <- max(c(infoMatch$vars[infoMatch[,varGroup]==pairGroups[1]]/infoMatch$vars[infoMatch[,varGroup]==pairGroups[2]],
                     infoMatch$vars[infoMatch[,varGroup]==pairGroups[2]]/infoMatch$vars[infoMatch[,varGroup]==pairGroups[1]]))

  return(list(stdzDiff = stdzDiff,
              ratioVars = ratioVars))

}

#' Balance indicators for binary variables
#' @keywords internal
balanceBinVar <- function(data, varBalance, match_id, varGroup, pairGroups, varWeights){


  varGroup_symbol <- as.symbol(varGroup)
  varWeights_symbol <- as.symbol(varWeights)

  #The "1" in the binary variables is considered the second level
  data$varBalanceNum <- (data[,varBalance] == (levels(data[,varBalance])[2]))*1

  #For difference in means, matched data
  dataMatch <- data[data[,varGroup] %in% pairGroups & !is.na(match_id), ]
  infoMatch <- dataMatch %>%
    dplyr::group_by(!!varGroup_symbol) %>%
    dplyr::summarize(means = stats::weighted.mean(x = .data$varBalanceNum, w = !!varWeights_symbol))

  #For variances, unmatched data
  dataUnm <- data[data[,varGroup] %in% pairGroups, ]
  infoUnm <- dataUnm %>%
    dplyr::group_by(!!varGroup_symbol) %>%
    dplyr::summarize(means = stats::weighted.mean(x = .data$varBalanceNum, w = !!varWeights_symbol)) %>%
    dplyr::mutate(vars = .data$means * (1 - .data$means))

  stdzDiff <- (infoMatch$means[infoMatch[,varGroup]==pairGroups[1]] - infoMatch$means[infoMatch[,varGroup]==pairGroups[2]])/(
               sqrt((infoUnm$vars[infoUnm[,varGroup]==pairGroups[1]] + infoUnm$vars[infoUnm[,varGroup]==pairGroups[2]])/2))

  return(list(stdzDiff=stdzDiff))

}

#' Balance indicators for categorical variables (more than 2 levels)
#' @keywords internal
balanceCatVar <- function(data, varBalance, match_id, varGroup, pairGroups, varWeights){

  varGroup_symbol <- as.symbol(varGroup)
  varBalance_symbol <- as.symbol(varBalance)
  varWeights_symbol <- as.symbol(varWeights)

  #For difference in means, matched data
  dataMatch <- data[data[,varGroup] %in% pairGroups & !is.na(match_id),]
  #Generate dummy variables
  dummyMatch <- as.data.frame(stats::model.matrix(stats::as.formula(paste("~ -1 + ", varBalance)),
                                           data = dataMatch))
  dataMatch <- cbind(dataMatch, dummyMatch)
  infoMatch <- dataMatch %>%
    dplyr::group_by(!!varGroup_symbol) %>%
    dplyr::summarize(dplyr::across(names(dummyMatch),
                                   ~stats::weighted.mean(x = .x, w = !!varWeights_symbol)),
                     .groups = "drop")


  p1 <- as.matrix(infoMatch[infoMatch[,varGroup] == pairGroups[1], names(dummyMatch)])
  p2 <- as.matrix(infoMatch[infoMatch[,varGroup] == pairGroups[2], names(dummyMatch)])

  p1 <- p1[1:(length(p1)-1)]
  p2 <- p2[1:(length(p2)-1)]

  #For variances, unmatched data
  dataUnm <- data[data[,varGroup] %in% pairGroups,]
  #Generate dummy variables
  dummyUnm <- as.data.frame(stats::model.matrix(stats::as.formula(paste("~ -1 + ", varBalance)),
                                           data = dataUnm))
  dataUnm <- cbind(dataUnm, dummyUnm)
  infoUnm <- dataUnm %>%
    dplyr::group_by(!!varGroup_symbol) %>%
    dplyr::summarize(dplyr::across(names(dummyUnm),
                                   ~stats::weighted.mean(x = .x, w = !!varWeights_symbol)),
                     .groups = "drop")


  p1unm <- as.matrix(infoUnm[infoUnm[,varGroup] == pairGroups[1], names(dummyUnm)])
  p2unm <- as.matrix(infoUnm[infoUnm[,varGroup] == pairGroups[2], names(dummyUnm)])

  p1unm <- p1unm[1:(length(p1unm)-1)]
  p2unm <- p2unm[1:(length(p2unm)-1)]


  Stemp <- (as.matrix(p1unm) %*% t(as.matrix(p1unm)) + as.matrix(p1unm) %*% t(as.matrix(p1unm)))/2
  diag(Stemp) <- 0
  S <- Stemp + (diag(p1unm*(1-p1unm)) + diag(p2unm*(1-p2unm)))/2

  stdzDiff <- sqrt(t(as.matrix(p1-p2)) %*% solve(S) %*% as.matrix(p1-p2))

  return(list(stdzDiff=stdzDiff))

}
