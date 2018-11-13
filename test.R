library(polymatching)

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


#Check balance (1) - different type of variables
#-----------------------------------------------
set.seed(123456)
dat <- generateData(c(100,200,200),
                    par = list(means = c(0,0,.5),
                               sds = c(1,3,1)))
dat$var1 <- factor(apply(rmultinom(nrow(dat), size = 1 , prob = c(1/3,1/3,1/3)), FUN = function(x){which(x==1)},2))
dat$var2 <- factor(rbinom(nrow(dat), size = 1, prob =.3))
dat$var3 <- rpois(nrow(dat),lambda = 100)
dat$var4 <- runif(nrow(dat))
dat$group <- factor(dat$group)

#Match (only on one variable)
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                     distance = "euclidean",
                     start = "1-2-3",
                     iterate = T, niter_max = 50, verbose = T)
dat$match_id <- result$match_id


#Plot distribution of matching variable
library(ggplot2)
ggplot(data = dat) + geom_density(aes(x=variable,
                                      colour = group,
                                      group = group))

ggplot(data = dat[!is.na(dat$match_id),]) +
  geom_density(aes(x=variable,
                   colour = group,
                   group = group))

resultBalance <- balance(group ~ variable + var1 + var2 + var3 + var4,
                        data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)


#Check balance (2) - several matching variables
#-----------------------------------------------
set.seed(123456)
dat <- data.frame(group = c(rep(1,100),rep(2,200),rep(3,200),rep(4,200)),
                  var1 = c(rnorm(100,0,1),rnorm(200,1,2),rnorm(200,-1,2),rnorm(200,0,3)),
                  var2 = c(rnorm(100,5,1),rnorm(200,6,2),rnorm(200,4,2),rnorm(200,5,3)),
                  var3 = c(rnorm(100,-5,1),rnorm(200,-6,2),rnorm(200,-4,2),rnorm(200,-5,3)))

result <- polymatch(formulaMatch = group ~ var1+var2+var3, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4",
                    iterate = T, niter_max = 50, verbose = T)

resultBalance <- balance(group ~  var1 + var2 + var3,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)

#Check balance (3) - our matched dataset
#-----------------------------------------------
load("C:/Users/natt03/Desktop/triplet matching/analysis/bestResult_3wayconstr_allData_1-2.Rdata")
dat <- result3wayConstr$matchedData

dat$income_1 <- factor(dat$income_1)
dat$income_2 <- factor(dat$income_2)
dat$income_3 <- factor(dat$income_3)
dat$income_4 <- factor(dat$income_4)
dat$pay1_1 <- factor(dat$pay1_1)
dat$pay1_2 <- factor(dat$pay1_2)
dat$pay1_3 <- factor(dat$pay1_3)
dat$pay1_4 <- factor(dat$pay1_4)
dat$pay1_5 <- factor(dat$pay1_5)
dat$pay1_6 <- factor(dat$pay1_6)
dat$nchs_1 <- factor(dat$nchs_1)
dat$nchs_2 <- factor(dat$nchs_2)
dat$nchs_3 <- factor(dat$nchs_3)
dat$nchs_4 <- factor(dat$nchs_4)
dat$nchs_5 <- factor(dat$nchs_5)
dat$nchs_6 <- factor(dat$nchs_6)
dat$chronic <- factor(dat$chronic)
dat$multiple_injury <- factor(dat$multiple_injury)
dat$FEMALE <- factor(dat$FEMALE)

#E.g. 1 - very long plot
resultBalance <- balance(HOSP_TRAUMA ~ AGE + iss +
                           income_1 + income_2 + income_3 + income_4 +
                           pay1_1 + pay1_2 + pay1_3 + pay1_4 + pay1_5 + pay1_6 +
                           nchs_1 + nchs_2 + nchs_3 + nchs_4 + nchs_5 + nchs_6 +
                           chronic + multiple_injury +
                           FEMALE,
                         data = dat, match_id = dat$indexMatch)
resultBalance

resultPlot <- plotBalance(resultBalance)
pdf("C:/Users/natt03/Desktop/triplet matching/analysis/balancePlot.pdf", height = 15, width = 6)
print(resultPlot[[1]])
dev.off()

#E.g. 2 - split balance in two
resultBalance1 <- balance(HOSP_TRAUMA ~ AGE + iss +
                           income_1 + income_2 + income_3 + income_4 +
                           pay1_1 + pay1_2 + pay1_3 + pay1_4 + pay1_5 + pay1_6,
                         data = dat, match_id = dat$indexMatch)
resultBalance2 <- balance(HOSP_TRAUMA ~ nchs_1 + nchs_2 + nchs_3 + nchs_4 + nchs_5 + nchs_6 +
                            chronic + multiple_injury + FEMALE,
                         data = dat, match_id = dat$indexMatch)
resultPlot1 <- plotBalance(resultBalance1)
resultPlot2 <- plotBalance(resultBalance2)

library(gridExtra)
library(ggplot2)
grid.arrange(resultPlot1[[1]] + labs(title=""),resultPlot2[[1]] + labs(title=""), ncol = 2)

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


