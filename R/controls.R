
#Check that type of inputs is appropriate
checkInputs <- function(formulaMatch, data, distance, start, iterate, niter_max, verbose){

  if( ! "formula" %in% class(formulaMatch) ) {
    stop("'formulaMatch' must be of class 'formula'")
  }

}

#Check coherence of inputs with data
checkData <- function(formulaMatch, data, start){

  varsFromFormula <- all.vars(formulaMatch)
  varGroup <- varsFromFormula[1]
  varsMatch <- varsFromFormula[2:length(varsFromFormula)]

  if( ! all(varsFromFormula %in% names(data)) ) {
    stop("The dataset provided must contain all the variables in 'formulaMatch'")
  }

  if("character" %in% class(start) & length(start)==1) {

    tabGroup <- table(data[,varGroup])

    if(start == "small.to.large") {

      vectorScheme <- names(sort(tabGroup))

    } else {

      vectorScheme <- unlist(strsplit(start, split = "-"))

      if(length(vectorScheme) != length(tabGroup)) {
        stop(paste0("'start' must contain the labels of ",length(tabGroup),
                    " groups separated by '-' (e.g.: '",
                    paste(names(sort(tabGroup)), collapse = "-"),"')"))
      }

      if(!all(sort(names(tabGroup)) == sort(vectorScheme))) {
        stop("The group labels in 'start' do not coincide with the group labels in the data")
      }

    }

  } else {

    vectorScheme <- NULL

  }

  if(length(start) > 1) {


    # Check that the starting point has 1 subject per group


  }

  return(list(varGroup = varGroup,
              varsMatch = varsMatch,
              vectorScheme = vectorScheme))
}
