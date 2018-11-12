#Check type of variable
typeVariable <- function(variable) {

  if(class(variable) %in% c("numeric","integer")) {

    type <- "continuous"

  } else {

    if(class(variable)=="factor") {

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

#Balance indicators for continuous variables
balanceContVar <- function(data, varBalance, match_id, varGroup, pairGroups){

  #For difference in means, matched data
  dataMeans <- data[data[,varGroup] %in% pairGroups & !is.na(match_id), ]
  means <- tapply(dataMeans[,varBalance], INDEX = dataMeans[,varGroup], FUN = mean)
  vars <- tapply(dataMeans[,varBalance], INDEX = dataMeans[,varGroup], FUN = stats::var)

  #For variances, unmatched data
  dataUnm <- data[data[,varGroup] %in% pairGroups, ]
  varsUnm <- tapply(dataUnm[,varBalance], INDEX = dataUnm[,varGroup], FUN = stats::var)

  stdzDiff <- (means[names(means)==pairGroups[1]] - means[names(means)==pairGroups[2]])/(
                sqrt((varsUnm[names(varsUnm)==pairGroups[1]] + varsUnm[names(varsUnm)==pairGroups[2]])/2))

  #Return the ratio of the variances >= 1
  ratioVars <- max(c(vars[names(vars)==pairGroups[1]]/vars[names(vars)==pairGroups[2]],
                     vars[names(vars)==pairGroups[2]]/vars[names(vars)==pairGroups[1]]))

  return(list(stdzDiff = stdzDiff,
              ratioVars = ratioVars))

}

#Balance indicators for binary variables
balanceBinVar <- function(data, varBalance, match_id, varGroup, pairGroups){

  #The "1" in the binary variables is considered the second level
  data$varBalanceNum <- (data[,varBalance] == (levels(data[,varBalance])[2]))*1

  #For difference in means, matched data
  dataMeans <- data[data[,varGroup] %in% pairGroups & match_id, ]
  means <- tapply(dataMeans$varBalanceNum, INDEX = dataMeans[,varGroup], FUN = mean)

  #For variances, unmatched data
  dataMeansUnm <- data[data[,varGroup] %in% pairGroups, ]
  meansUnm <- tapply(dataMeansUnm$varBalanceNum, INDEX = dataMeansUnm[,varGroup], FUN = mean)
  vars <- meansUnm*(1-meansUnm)

  stdzDiff <- (means[names(means)==pairGroups[1]] - means[names(means)==pairGroups[2]])/(
               sqrt((vars[names(vars)==pairGroups[1]] + vars[names(vars)==pairGroups[2]])/2))

  return(list(stdzDiff=stdzDiff))

}

#Balance indicators for categorical variables (more than 2 levels)
balanceCatVar <- function(data, varBalance, match_id, varGroup, pairGroups){

  #For difference in means, matched data
  dataMeans <- data[data[,varGroup] %in% pairGroups & !is.na(match_id), ]

  p1 <- prop.table(table(dataMeans$var1[dataMeans$group %in% pairGroups[1]]))
  p2 <- prop.table(table(dataMeans$var1[dataMeans$group %in% pairGroups[2]]))

  p1 <- p1[1:(length(p1)-1)]
  p2 <- p2[1:(length(p2)-1)]

  #For variances, unmatched data
  p1unm <- prop.table(table(dataMeans$var1[dataMeans$group %in% pairGroups[1]]))
  p2unm <- prop.table(table(dataMeans$var1[dataMeans$group %in% pairGroups[2]]))

  p1unm <- p1unm[1:(length(p1unm)-1)]
  p2unm <- p2unm[1:(length(p2unm)-1)]

  Stemp <- (as.matrix(p1unm) %*% t(as.matrix(p1unm)) + as.matrix(p1unm) %*% t(as.matrix(p1unm)))/2
  diag(Stemp) <- 0
  S <- Stemp + (diag(p1unm*(1-p1unm)) + diag(p2unm*(1-p2unm)))/2

  stdzDiff <- sqrt(t(as.matrix(p1-p2)) %*% solve(S) %*% as.matrix(p1-p2))

  return(list(stdzDiff=stdzDiff))

}
