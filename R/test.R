source("R/polymatch.R")
source("R/controls.R")
source("R/matching_functions.R")
source("R/balance.R")
source("R/balance_functions.R")

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


#Check iterations
#----------------
set.seed(123456)
dat <- generateData(c(100,100,100,100,100,100))

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                     distance = "euclidean",
                     start = "1-2-3-4-5-6",
                     iterate = T, niter_max = 50, verbose = T)

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4-5-6",
                    iterate = F, niter_max = 50, verbose = T)

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4-5-6",
                    iterate = F, niter_max = 3, verbose = T)

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4-5-6",
                    iterate = T, niter_max = 3, verbose = T)

#Verbose
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4-5-6",
                    iterate = T, niter_max = 50, verbose = F)

#Check distances
#----------------
set.seed(123456)
dat <- data.frame(group = c(rep(1,10), rep(2,50), rep(3,50)),
                  var1 = rnorm(110),
                  var2 = rnorm(110)*1000)

resultE <- polymatch(formulaMatch = group ~ var1 + var2, data = dat,
                    distance = "euclidean",
                    start = "1-2-3",
                    iterate = T, niter_max = 50, verbose = T)

resultM <- polymatch(formulaMatch = group ~ var1 + var2, data = dat,
                     distance = "mahalanobis",
                     start = "1-2-3",
                     iterate = T, niter_max = 50, verbose = T)

dat$match_idE <- resultE$match_id
dat$match_idM <- resultM$match_id

dat[dat$match_idE %in% 1,]
dat[dat$match_idM %in% 1,]


#Check balance
#----------------
set.seed(123456)
dat <- generateData(c(100,200,200),
                    par = list(means = c(0,0,2),
                               sds = c(1,3,1)))
dat$var1 <- factor(apply(rmultinom(nrow(data), size = 1 , prob = c(1/3,1/3,1/3)), FUN = function(x){which(x==1)},2))
dat$var2 <- factor(rbinom(nrow(data), size = 1, prob =.3))
dat$var3 <- rpois(nrow(data),lambda = 100)
dat$var4 <- runif(nrow(data))


result <- polymatch(formulaMatch = group ~ variable, data = dat,
                     distance = "euclidean",
                     start = "1-2-3",
                     iterate = T, niter_max = 50, verbose = T)

resultBalance <- balance(group ~ variable + var1 + var2 + var3 + var4,
                        data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plot.balanceCondOptMatch(resultBalance)


#Example with some numerical weird issues
#----------------------------------------

#When there is a very large overlap among the observations, the
#optimal solution identified by package 'optmatch' is NOT always the same.
#If the labels of the groups are held constant, the solution seems to be always the same.
#However, if I relabel the name of the groups (e.g., "1-2-3-4-5" as "a-b-c-d-e") I might find
#something slightly different. This bizarre result is reproduced below.

#I give different names to the groups and run our algorithm. The optimal solution "seem" to be the same
#in terms of total distance (it's probably a result of numerical rounding). BUT, the number of iterations is different,
#even though we have the same starting point. AND, the final matched set is not the same in the two cases!

#Run on R version 3.5.1
set.seed(123456)
dat <- generateData(c(100,100,100,100,100))

dat$groupN <- dat$group
dat$groupF <- factor(dat$group, levels = c(1:5),
                    labels = c("aaaa","asda","ACADDDDD","jkasdljads_dahskj","aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"))

resultF <- polymatch(formulaMatch = groupF ~ variable, data = dat,
                distance = "euclidean",
                start = "aaaa-asda-ACADDDDD-jkasdljads_dahskj-aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                iterate = T, niter_max = 50, verbose = T)

resultN <- polymatch(formulaMatch = groupN ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4-5",
                    iterate = T, niter_max = 50, verbose = T)

#Different number of iterations
resultF$niter
resultN$niter

#Different matched sets (I expect only 0s and 5s if the matched sets would be the same)
table(table(resultF$match_id, resultN$match_id))

#Example where anchor matching is poor (Bo Lu)
#---------------------------------------------

dat <- data.frame(group = c("A","A","B","B","B","C","C","C"),
                  variable = c(5,10,1,8,11,2,8,13))

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean", start = "small.to.large",
                    iterate = T, niter_max = 50, verbose = T)

dat$match_id <- result$match_id

#Optimal solution identified (5-8-8 and 10-11-13)
dat


