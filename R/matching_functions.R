#'Distance between two units
#' @keywords internal
pairwiseDistance <- function(A, B, distance, Sigma) {

  if(distance=="euclidean") {
    output <- sqrt(sum((A-B)^2))
  }

  if(distance=="mahalanobis") {
    output <- sqrt(stats::mahalanobis(x = A, center = B, cov = Sigma))
  }

  return(output)
}

#'Given a set of units, compute all the pairwise distances.
#' @keywords internal
computePairwiseDistances <- function(dataPolygon,  distance, Sigma) {

  dataPolygon <- as.matrix(dataPolygon)

  pairGroups <- utils::combn(1:nrow(dataPolygon), 2)

  sumPairwDist <- 0

  for(indexPair in 1:ncol(pairGroups)) {
    pairGroupsTemp <- pairGroups[,indexPair]
    sumPairwDist <- (sumPairwDist +
                     pairwiseDistance(dataPolygon[pairGroupsTemp[1],],
                                      dataPolygon[pairGroupsTemp[2],],
                                      distance = distance,
                                      Sigma = Sigma))
  }

  return(sumPairwDist)
}

#' Compute total distance of a matched sample
#' @keywords internal
evaluateMatching <- function(data, varIndexMatch, varsMatch, distance, Sigma) {

  withinMatchDistances <- by(data[,varsMatch],
                             INDICES = data[,varIndexMatch],
                             FUN = computePairwiseDistances,
                             distance = distance,
                             Sigma = Sigma)

  dfTempToExport <- data.frame(indexMatch = names(withinMatchDistances),
                               distance = as.vector(withinMatchDistances),
                               stringsAsFactors = F)

  names(dfTempToExport)[names(dfTempToExport) == "indexMatch"] <- varIndexMatch

  return(list(total_distance = sum(withinMatchDistances),
              distanceByMatch = dfTempToExport))
}


#' Function implementing conditional optimal matching.
#'
#' Optimally matches two sets of IDs of matched units--varIndexMatch1 and varIndexMatch2--with respect
#' to matching variables varsMatch.
#'
#' @keywords internal
condOptMatching <- function(data, varIndexMatch1, varIndexMatch2,
                            varsMatch, varGroup,
                            distance, Sigma) {


  #Local function for conditional matching:
  #----------------------------------------

  applyPersonalDistance <- function(index, data, z) {

    #browser()

    indexDf <- as.data.frame(index, stringsAsFactors = F)
    groupTreated <- unique(data[z,varGroup])
    groupControls <- unique(data[!z,varGroup])
    names(indexDf) <- paste("group",c(groupTreated,groupControls), sep ="")

    #I can't pass additional arguments to the function 'applyPersonalDistance' through 'match_on',
    # 'match_on' passes ... to internal functions. To access to other objects, I need to take them
    # directly from the parent envirormnet, which is the environment where the function 'applyPersonalDistance' lives.
    #Grab objects from "outside the function":
    envDataAll <- environment(applyPersonalDistance)
    dataAll <- get("dataAll",envir = envDataAll)
    Sigma <- get("Sigma",envir = envDataAll)
    distance <- get("distance",envir = envDataAll)

    #Wide format for:

    # - group 1
    longIndexGroup1 <- dataAll[!is.na(dataAll[,varIndexMatch1]), c(varGroup, varIndexMatch1)]
    longIndexGroup1$value <- row.names(longIndexGroup1)
    wideIndexGroup1 <- stats::reshape(longIndexGroup1,
                               direction = "wide",
                               idvar = varIndexMatch1,
                               timevar = varGroup)
    wideIndexGroup1[,varIndexMatch1] <- NULL
    names(wideIndexGroup1) <- gsub("value.","group",names(wideIndexGroup1))

    # - group 2
    longIndexGroup2 <- dataAll[!is.na(dataAll[,varIndexMatch2]), c(varGroup, varIndexMatch2)]
    longIndexGroup2$value <- row.names(longIndexGroup2)
    wideIndexGroup2 <- stats::reshape(longIndexGroup2,
                               direction = "wide",
                               idvar = varIndexMatch2,
                               timevar = varGroup)
    wideIndexGroup2[,varIndexMatch2] <- NULL
    names(wideIndexGroup2) <- gsub("value.","group",names(wideIndexGroup2))

    #Merge the groups not imputed in the function to indexDf:

    # - group 1
    if(length(groups1)>=2) {
      matchingVector1 <- match(indexDf[, paste("group",groups1[1],sep="")],
                               wideIndexGroup1[, paste("group", groups1[1], sep="")])
      indexDf[, paste("group",
                      groups1[2:length(groups1)], sep = "") ] <- wideIndexGroup1[matchingVector1,
                                                                                 paste("group", groups1[2:length(groups1)], sep="")]
    }
    # - group 2
    if(length(groups2)>=2) {
      matchingVector2 <- match(indexDf[, paste("group",groups2[1],sep="")],
                               wideIndexGroup2[, paste("group", groups2[1], sep="")])
      indexDf[, paste("group",
                      groups2[2:length(groups2)], sep = "") ] <- wideIndexGroup2[matchingVector2,
                                                                                 paste("group", groups2[2:length(groups2)], sep="")]
    }

    #Merge with matching variables for each group
    varsByGroup <- paste("group",outer(c(groups1, groups2), varsMatch, FUN = paste, sep="_"),sep="")

    indexDf[,varsByGroup] <- NA
    for(groupTemp in c(groups1, groups2)) {
      indexDf[, paste("group",
                      groupTemp, "_",
                      varsMatch, sep = "")] <- dataAll[indexDf[,paste("group",
                                                                      groupTemp,
                                                                      sep = "")],varsMatch]
    }

    #Sum the pairwise distances within matched sets
    pairGroups <- utils::combn(c(groups1, groups2), 2)
    distances <- rep(0, nrow(indexDf))

    #For Mahalanobis distance, need the inverse of Sigma
    if(distance == "mahalanobis") {
      SigmaInv <- try(chol2inv(chol(Sigma)), silent = T)
      if(class(SigmaInv)=="try-error") {
        stop("Problems in the computation of Mahalanobis distance, error inverting vcov matrix. Try using distance='euclidean'. ")
      }
    }

    for(indexPair in 1:ncol(pairGroups)) {

      pairGroupsTemp <- pairGroups[,indexPair]

      matrixTemp <- as.matrix((indexDf[,paste("group",pairGroupsTemp[1],"_",varsMatch,sep="")] -
                                 indexDf[,paste("group",pairGroupsTemp[2],"_",varsMatch,sep="")]))

      if(distance == "euclidean") {
        distances <- (distances + sqrt(rowSums(matrixTemp^2)))
      }
      if(distance == "mahalanobis") {
        distances <- (distances + sqrt(rowSums((matrixTemp %*% SigmaInv) * matrixTemp)))
      }

    }

    return(distances)
  }

  ############################################################################################

  groups1 <- sort(unique(data[!is.na(data[,varIndexMatch1]), varGroup]))
  groups2 <- sort(unique(data[!is.na(data[,varIndexMatch2]), varGroup]))

  #In the selection, selecting groups1[1] and groups2[1] is arbitrary.
  #However, the same choice made here must be reported also into the function applyPersonalDistance
  selectionToMatch <- ( (data[,varGroup] %in% groups1[1] & !is.na(data[,varIndexMatch1])) |
                          (data[,varGroup] %in% groups2[1] & !is.na(data[,varIndexMatch2])) )

  #Generate new variable that define the "bipartite" groups
  data$groupNew <- NA
  data$groupNew[data[,varGroup] %in% groups1] <- paste(groups1, collapse = "")
  data$groupNew[data[,varGroup] %in% groups2] <- paste(groups2, collapse = "")
  data$groupNew <- factor(data$groupNew,
                          levels = c(paste(groups1, collapse = ""),
                                     paste(groups2, collapse = "")))
  tabTemp <- table(data$groupNew[selectionToMatch])
  data$groupNew <- stats::relevel(data$groupNew,
                                  names(tabTemp)[which.max(tabTemp)])

  #Make treatment variable binary. Package optmatch deprecated factor type for treatment variables.
  data$groupNew <- (data$groupNew == (levels(data$groupNew)[2]))*1

  #Global variables to be used in the function 'applyPersonalDistance'
  dataAll <-  data

  #browser()
  resultDistance <- optmatch::match_on(applyPersonalDistance,
                             z = data$groupNew[selectionToMatch],
                             data = data[selectionToMatch,])
  #optmatch::caliper(resultDistance, width = 0.649, values = TRUE)
  #The possible caliper is a direct truncation of all the distances within that value

  resultMatch <- optmatch::pairmatch(resultDistance,
                           controls = 1,
                           data = data[selectionToMatch,])

  data$indexResultMatch <- NA
  data$indexResultMatch[selectionToMatch] <- resultMatch

  # Match all groups
  #-----------------

  indexMatchesInGroup1 <- !is.na(data$indexResultMatch) & !is.na(data[,varIndexMatch1])

  data$indexMatch <- NA
  matchedCounter <- 0

  for (i in which(indexMatchesInGroup1)) {

    matchedCounter <- matchedCounter + 1
    idGroup1 <- data[i,varIndexMatch1]
    idGroup2 <- data[data$indexResultMatch %in% data$indexResultMatch[i] &
                       !is.na(data[,varIndexMatch2]), varIndexMatch2]

    data$indexMatch[data[,varIndexMatch1] %in% idGroup1 |
                      data[,varIndexMatch2] %in% idGroup2] <- matchedCounter

  }

  resultEvaluation <- evaluateMatching(data,
                                       varIndexMatch = "indexMatch",
                                       varsMatch = varsMatch,
                                       distance = distance,
                                       Sigma = Sigma)

  return(list(total_distance = resultEvaluation$total_distance,
              match_id = data$indexMatch))
}
