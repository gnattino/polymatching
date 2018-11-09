source("R/polymatching.R")
source("R/controls.R")
source("R/matching_functions.R")

##########################################
# Function to generate data to polymatch #
##########################################

generateData <- function(nVect, par = NULL) {

  nGroups <- length(nVect)

  #If null: all data generated from N(0,1)
  if(is.null(par)) {
    par <- list()
    par$means <- rep(0,nGroups)
    par$sds <- rep(1,nGroups)
  }


  result <- data.frame(group = unlist(mapply(1:length(nVect),
                                             each = nVect,
                                             FUN = rep, SIMPLIFY = F)),
                       variable = NA)

  for(group in 1:nGroups) {

    result[result$group %in% group, "variable"] <- rnorm(nVect[group],
                                                         mean = par$means[group],
                                                         sd = par$sd[group])

  }

  return(result)

}

#########
# Check #
#########

dat <- generateData(c(100,100,100,100,100))

dat$groupN <- dat$group
dat$groupF <- factor(dat$group, levels = c(1:5),
                    labels = c("aaaa","asda","ACADDDDD","jkasdljads_dahskj","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"))

data <- dat
dataBKP <- dat


dat$variable <- dat$variable*10000000

result <- polymatch(formulaMatch = groupF ~ variable, data = dat,
                distance = "euclidean", start = "aaaa-asda-ACADDDDD-jkasdljads_dahskj-aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                iterate = T, niter_max = 50, verbose = T)

result <- polymatch(formulaMatch = groupN ~ variable, data = dat,
                    distance = "euclidean", start = "1-2-3-4-5",
                    iterate = T, niter_max = 50, verbose = T)


dataCondOpt$groupF <- factor(dataCondOpt$groupN, levels = c(1:5),
                     labels = c("aaaa","asda","ACADDDDD","jkasdljads_dahskj","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"))
dataCondOpt$groupF <- as.character(dataCondOpt$groupF)


dataCondOpt$groupN <- factor(dataCondOpt$groupF, labels = c(1:5),
                             levels = c("aaaa","asda","ACADDDDD","jkasdljads_dahskj","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"))
dataCondOpt$groupN <- as.numeric(as.character(dataCondOpt$groupN))


resultIter <- condOptMatching(data = dataCondOpt,
                              varIndexMatch1 = "indexMatchIter1",
                              varIndexMatch2 = "indexMatchIter2",
                              varsMatch = "variable",
                              varGroup = "groupF",
                              distance = "euclidean",
                              Sigma = NULL,
                              niter = 1000)

sprintf("%.3f",resultIter$total_distance)

resultIter <- condOptMatching(data = dataCondOpt,
                              varIndexMatch1 = "indexMatchIter1",
                              varIndexMatch2 = "indexMatchIter2",
                              varsMatch = "variable",
                              varGroup = "groupN",
                              distance = "euclidean",
                              Sigma = NULL,
                              niter = 1000)


#Example from dr. Lu
dat <- data.frame(#group = c("A","A","B","B","B","C","C","C"),
                  group = c(1,1,2,2,2,3,3,3),
                  variable = c(5,10,1,8,11,2,8,13))

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean", start = "small.to.large",
                    iterate = T, niter_max = 50, verbose = T)

#Understand:
# - why different rseults when change labels of names?
# - why can't find optimal solution in example?

