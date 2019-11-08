
#' Check that type of inputs is appropriate
#' @keywords internal
checkInputs <- function(formulaMatch, start, data, distance, exactMatch, iterate, niter_max, verbose){

  if( ! "formula" %in% class(formulaMatch) ) {
    stop("'formulaMatch' must be of class 'formula'")
  }

  if( !is.null(exactMatch)) {
    if( ! "formula" %in% class(exactMatch) ) {
      stop("'exactMatch' must be of class 'formula'")
    }
  }

  if( length(start)>1 & length(start)<nrow(data)) {
    stop("'start' must be either a character string or a vector of length nrow(data)")
  }

}

#' Check coherence of inputs with data
#' @keywords internal
checkData <- function(formulaMatch, start, data, exactMatch, checkOnePerGroup = TRUE){
  #browser()

  varsFromFormula <- all.vars(formulaMatch)
  varGroup <- varsFromFormula[1]
  varsMatch <- varsFromFormula[2:length(varsFromFormula)]

  tabGroup <- table(data[,varGroup])

  if( ! all(varsFromFormula %in% names(data)) ) {
    stop("The dataset provided must contain all the variables in 'formulaMatch'")
  }

  if( sum(is.na(data[,varGroup])) > 0 ) {
    stop(paste0("The variable with the treatment groups (",varGroup,") cannot have missing values"))
  }

  if("character" %in% class(start) & length(start)==1) {

    if(start == "small.to.large") {

      vectorSchemeStart <- names(sort(tabGroup))

    } else {

      vectorSchemeStart <- unlist(strsplit(start, split = "-"))

      if(length(vectorSchemeStart) != length(tabGroup)) {
        stop(paste0("'start' must contain the labels of ",length(tabGroup),
                    " groups separated by '-' (e.g.: '",
                    paste(names(sort(tabGroup)), collapse = "-"),"')"))
      }

      if(!all(sort(names(tabGroup)) == sort(vectorSchemeStart))) {
        stop("The group labels in 'start' do not coincide with the group labels in the data")
      }

    }

  } else {

    vectorSchemeStart <- NULL

  }

  if(length(start) > 1) {
    
    if(checkOnePerGroup == TRUE) {
      
      # Check that the starting point has exactly 1 subject per group
      tabStartGroup <- table(start, data[,varGroup])
      if(!all(tabStartGroup==1)) {
        stop("The matched sets provided in 'start' must have exactly one subject per group")
      }
      
    }

    
    #Check that all the subjects in the smallest group(s) are matched
    smallestGroups <- names(tabGroup)[tabGroup==min(tabGroup)]
    if( sum(is.na(start[data[,varGroup] %in% smallestGroups]))>0) {
      warning("In the matched sample provided, some units of the smallest group(s) are NOT matched")
    }


  }

  if(!is.null(exactMatch)) {

    varsExactMatch <- all.vars(exactMatch)

    #Check: all the variables are in the dataset
    if( ! all(varsExactMatch %in% names(data)) ) {
      stop("The dataset provided must contain all the variables in 'exactMatch'")
    }

    #Check: if there is an initial dataset, it must be exactly matched on the variables.
    if(length(start) > 1) {

      for(varExactMatchTemp in varsExactMatch) {

        if(!all(rowSums(table(start,data[,varExactMatchTemp])==length(tabGroup))==1)) {
          stop(paste0("The matched sample provided is not exactly matched on '",varExactMatchTemp,"'"))
        }

      }

    }



  } else {
    varsExactMatch <- NULL
  }

  return(list(varGroup = varGroup,
              varsMatch = varsMatch,
              vectorSchemeStart = vectorSchemeStart,
              varsExactMatch = varsExactMatch))
}
