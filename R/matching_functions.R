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
computePairwiseDistances <- function(data,  distance, Sigma, varGroup, varsMatch, dat_stdzDistances) {

  groups <- dplyr::pull(data, varGroup)
  
  dataPolygon <- as.matrix(data[, varsMatch])

  pairGroups <- utils::combn(1:nrow(dataPolygon), 2)

  sumPairwDist <- 0

  for(indexPair in 1:ncol(pairGroups)) {
    
    pairGroupsTemp <- pairGroups[,indexPair]
    
    sortedGroups <- sort(c(groups[pairGroupsTemp[1]], groups[pairGroupsTemp[2]]))
    
    factor_stdzDistances <- dat_stdzDistances$factor[ dat_stdzDistances$group_1 == sortedGroups[1] &
                                                        dat_stdzDistances$group_2 == sortedGroups[2] ]
    
    
    sumPairwDist <- (sumPairwDist +
                     pairwiseDistance(dataPolygon[pairGroupsTemp[1],],
                                      dataPolygon[pairGroupsTemp[2],],
                                      distance = distance,
                                      Sigma = Sigma)/factor_stdzDistances)
  }

  return(sumPairwDist)
}

#' Compute total distance of a matched sample
#' @keywords internal
evaluateMatching <- function(data, varIndexMatch, varGroup, varsMatch, distance, Sigma, dat_stdzDistances) {

  withinMatchDistances <- by(data[,c(varGroup, varsMatch)],
                             INDICES = data[,varIndexMatch],
                             FUN = computePairwiseDistances,
                             distance = distance,
                             Sigma = Sigma,
                             varGroup = varGroup,
                             varsMatch = varsMatch,
                             dat_stdzDistances = dat_stdzDistances)

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
                            distance, Sigma,
                            varsExactMatch,
                            k,
                            dat_stdzDistances) {


  #Local function for conditional matching:
  #----------------------------------------

  applyPersonalDistance <- function(index, data, z) {

    #I can't pass additional arguments to the function 'applyPersonalDistance' through 'match_on',
    # 'match_on' passes ... to internal functions. To access to other objects, I need to take them
    # directly from the parent envirormnet, which is the environment where the function 'applyPersonalDistance' lives.
    #Grab objects from "outside the function":
    envExternalFun <- environment(applyPersonalDistance)
    dataAll <- get("dataAll",envir = envExternalFun)
    Sigma <- get("Sigma",envir = envExternalFun)
    distance <- get("distance",envir = envExternalFun)
    varIndexMatch1 <- get("varIndexMatch1",envir = envExternalFun)
    varIndexMatch2 <- get("varIndexMatch2",envir = envExternalFun)
    groups1 <- get("groups1",envir = envExternalFun)
    groups2 <- get("groups2",envir = envExternalFun)
    dat_stdzDistances <- get("dat_stdzDistances",envir = envExternalFun)

    #Generate indexDf is a matrix with all the treated-control pairs as rows.
    indexDf <- as.data.frame(index, stringsAsFactors = F)
    groupTreated <- unique(data$groupNewFunction[z])
    groupControls <- unique(data$groupNewFunction[!z])
    names(indexDf) <- paste("group",c(groupTreated,groupControls), sep ="")

    #Wide format for:

    # - group 1
    longIndexGroup1 <- dataAll[!is.na(dataAll[,varIndexMatch1]), c("groupNewFunction", varIndexMatch1)]
    longIndexGroup1$value <- row.names(longIndexGroup1)
    wideIndexGroup1 <- stats::reshape(longIndexGroup1,
                               direction = "wide",
                               idvar = varIndexMatch1,
                               timevar = "groupNewFunction")
    wideIndexGroup1[,varIndexMatch1] <- NULL
    names(wideIndexGroup1) <- gsub("value.","group",names(wideIndexGroup1))

    # - group 2
    longIndexGroup2 <- dataAll[!is.na(dataAll[,varIndexMatch2]), c("groupNewFunction", varIndexMatch2)]
    longIndexGroup2$value <- row.names(longIndexGroup2)
    wideIndexGroup2 <- stats::reshape(longIndexGroup2,
                               direction = "wide",
                               idvar = varIndexMatch2,
                               timevar = "groupNewFunction")
    wideIndexGroup2[,varIndexMatch2] <- NULL
    names(wideIndexGroup2) <- gsub("value.","group",names(wideIndexGroup2))

    #Merge the groups not imputed in the function to indexDf:

    # - group 1
    if(length(groups1)>=2) {
      matchingVector1 <- match(indexDf[, paste("group", groups1[1],sep="")],
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
      if("try-error" %in% class(SigmaInv) ) {
        stop("Problems in the computation of Mahalanobis distance, error inverting vcov matrix. Try using distance='euclidean'. ")
      }
    }

    for(indexPair in 1:ncol(pairGroups)) {

      pairGroupsTemp <- pairGroups[,indexPair]
      
      #Original name of groups
      originalNamesGroupsTemp <- sort(gsub("\\.[^.]*$", "", pairGroupsTemp))
      factor_stdzDistances <- dat_stdzDistances$factor[ dat_stdzDistances$group_1 == originalNamesGroupsTemp[1] &
                                                          dat_stdzDistances$group_2 == originalNamesGroupsTemp[2] ]
      
      matrixTemp <- as.matrix((indexDf[,paste("group",pairGroupsTemp[1],"_",varsMatch,sep="")] -
                                 indexDf[,paste("group",pairGroupsTemp[2],"_",varsMatch,sep="")]))

      if(distance == "euclidean") {
        distances <- (distances + sqrt(rowSums(matrixTemp^2))/factor_stdzDistances)
      }
      if(distance == "mahalanobis") {
        distances <- (distances + sqrt(rowSums((matrixTemp %*% SigmaInv) * matrixTemp))/factor_stdzDistances)
      }

    }

    return(distances)
  }

  ############################################################################################

  #Modify groups' name: when group X have k subjects matched,
  # define groups X.1, X.2, ..., X.k. Subjects are randomly allocated
  # to groups. This is not problematic as we compute all pairwise distances
  # so it does not matter which subject is in which group.
  data <- data %>%
    dplyr::group_by(!!as.symbol(varGroup), !!as.symbol(varIndexMatch1)) %>%
    dplyr::mutate(id_rep_group_1 = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!as.symbol(varGroup), !!as.symbol(varIndexMatch2)) %>%
    dplyr::mutate(id_rep_group_2 = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      groupNewFunction = dplyr::case_when(
        !is.na(!!as.symbol(varIndexMatch1)) ~ paste0(!!as.symbol(varGroup), ".", id_rep_group_1),
        !is.na(!!as.symbol(varIndexMatch2)) ~ paste0(!!as.symbol(varGroup), ".", id_rep_group_2)
      )
      ) %>%
    as.data.frame()
  
  data$id_rep_group_1 <- NULL
  data$id_rep_group_2 <- NULL
  
  groups1 <- sort(unique(data$groupNewFunction[!is.na(data[,varIndexMatch1])]))
  groups2 <- sort(unique(data$groupNewFunction[!is.na(data[,varIndexMatch2])]))

  #In the selection, selecting groups1[1] and groups2[1] is arbitrary.
  #However, the same choice made here must be reported also into the function applyPersonalDistance
  selectionToMatch <- ( (data$groupNewFunction %in% groups1[1] & !is.na(data[,varIndexMatch1])) |
                          (data$groupNewFunction %in% groups2[1] & !is.na(data[,varIndexMatch2])) )

  #Generate new variable that define the "bipartite" groups
  data$groupNew <- NA
  data$groupNew[data$groupNewFunction %in% groups1] <- paste(groups1, collapse = "")
  data$groupNew[data$groupNewFunction %in% groups2] <- paste(groups2, collapse = "")
  data$groupNew <- factor(data$groupNew,
                          levels = c(paste(groups1, collapse = ""),
                                     paste(groups2, collapse = "")))
  tabTemp <- table(data$groupNew[selectionToMatch])
  data$groupNew <- stats::relevel(data$groupNew,
                                  names(tabTemp)[which.max(tabTemp)])

  #Make treatment variable binary. Package optmatch deprecated factor type for treatment variables.
  data$groupNew <- (data$groupNew == (levels(data$groupNew)[2]))*1

  #Global variable to be used in the function 'applyPersonalDistance'
  dataAll <-  data

  #If there are some variables to match exactly on, do that:
  if(!is.null(varsExactMatch)) {

    formulaExactMatch <- stats::as.formula(paste0("groupNew ~ ",paste(varsExactMatch, collapse = "+")))

    resultExactMatch <- optmatch::exactMatch(formulaExactMatch, data = data[selectionToMatch,])

    resultDistance <- optmatch::match_on(applyPersonalDistance,
                                         z = data$groupNew[selectionToMatch],
                                         data = data[selectionToMatch,],
                                         within = resultExactMatch)

    } else {
      #Otherwise, matching without exact matching constraints:
      resultDistance <- optmatch::match_on(applyPersonalDistance,
                                           z = data$groupNew[selectionToMatch],
                                           data = data[selectionToMatch,])

  }

  #We can also think of adding the caliper.
  #optmatch::caliper(resultDistance, width = 0.649, values = TRUE)
  #The possible caliper is a direct truncation of all the distances within that value

  resultMatch <- optmatch::pairmatch(resultDistance,
                                     controls = k,
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
                                       varGroup = varGroup,
                                       varsMatch = varsMatch,
                                       distance = distance,
                                       Sigma = Sigma,
                                       dat_stdzDistances = dat_stdzDistances)

  return(list(total_distance = resultEvaluation$total_distance,
              match_id = data$indexMatch))
}

#'Compute factor to standardize distances when matching multiple subjects in one group
#' @keywords internal
stdzDistances <- function(vectorK, withinGroupDist) {
  
  groups <- names(vectorK)
  dat_stdzDistances <- NULL
  
  #Distances between groups
  #-------------------------
  pairs_groups <- utils::combn(groups, 2)
  
  for(j in 1:ncol(pairs_groups)) {
    
    sorted_groups <- sort(pairs_groups[,j])
    
    dat_stdzDistances_iter <- data.frame(
      group_1 = sorted_groups[1],
      group_2 = sorted_groups[2],
      factor = as.numeric(vectorK[sorted_groups[1]]) * as.numeric(vectorK[sorted_groups[2]]),
      stringsAsFactors = FALSE
      )
    
    dat_stdzDistances <- rbind(dat_stdzDistances,
                               dat_stdzDistances_iter)
    
  }
  
  #Distances within groups
  #------------------------
  for(i in 1:length(groups)) {
    
    dat_stdzDistances_iter <- data.frame(
      group_1 = groups[i],
      group_2 = groups[i],
      factor = ifelse(withinGroupDist, 
                      choose(as.numeric(vectorK[groups[i]]), 2),
                      Inf),
      stringsAsFactors = FALSE
    )
    
    dat_stdzDistances <- rbind(dat_stdzDistances,
                               dat_stdzDistances_iter)
    
  }
  
  return(dat_stdzDistances)
}
