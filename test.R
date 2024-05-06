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

#Check starting points
#---------------------
set.seed(123456)
dat <- generateData(c(80,50,90,75,200,100))

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "small.to.large",
                    iterate = T, niter_max = 50, verbose = T)

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "5-6-3-4-1-2",
                    iterate = T, niter_max = 50, verbose = T)

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4-5-6",
                    iterate = T, niter_max = 50, verbose = T)

dat$match_id <- NA
for(group in unique(dat$group)) {
  dat$match_id[dat$group %in% group][sample(sum(dat$group %in% group),50)] <- 1:50
}

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = dat$match_id,
                    iterate = T, niter_max = 50, verbose = T)

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


#Exact match (1) - everything goes fine
#----------------------------------------
set.seed(123456)
dat <- generateData(c(100,200,300,400),
                    par = list(means = c(0,0,0,.5),
                               sds = c(1,1,3,1)))
dat$var1 <- factor(apply(rmultinom(nrow(dat), size = 1 , prob = c(1/10,2/10,3/10,4/10)), FUN = function(x){which(x==1)},2))
dat$var2 <- factor(rbinom(nrow(dat), size = 1 , prob = c(1/10,9/10)), levels = c(0,1), labels = c("A","B"))

#I have one exact match for each of the subjects in the smallest group
table(dat$group,dat$var1, dat$var2, dnn = c("group", "var1","var2"))

#Match on one variable and exact match on other 2
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4",
                    exactMatch = ~var1+var2,
                    iterate = T, niter_max = 50, verbose = T)

resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)


#Match only on one variable, exact match only on one variable
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4",
                    exactMatch = ~var1,
                    iterate = T, niter_max = 50, verbose = T)
resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)

#Match only on one variable, no exact match
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4",
                    iterate = T, niter_max = 50, verbose = T)
resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)
#As expected, total distance is way smaller

#Bad choice of order
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "4-3-2-1",
                    exactMatch = ~var1+var2,
                    iterate = T, niter_max = 50, verbose = T)
resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)
#As expected, total distance is way smaller

#Exact match (2) - Exact match does not exists
#----------------------------------------
set.seed(123456)
dat <- generateData(c(100,100,100,100),
                    par = list(means = c(0,0,0,.5),
                               sds = c(1,1,3,1)))
dat$var1 <- factor(apply(rmultinom(nrow(dat), size = 1 , prob = c(1/10,2/10,3/10,4/10)), FUN = function(x){which(x==1)},2))
dat$var2 <- factor(rbinom(nrow(dat), size = 1 , prob = c(1/10,9/10)), levels = c(0,1), labels = c("A","B"))

#Now I dow't have one exact match for each of the subjects in the smallest group
table(dat$group,dat$var1, dat$var2, dnn = c("group", "var1","var2"))
#subjects that can be matched exactly are the minimum by column:
4+6+9+16+2+5+12+18

#If matching exactly only on var1:
table(dat$group,dat$var1, dnn = c("group", "var1"))
7+18+26+36


#If exact match does not exists, not all of the subjects are matched
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4",
                    exactMatch = ~var1+var2,
                    iterate = T, niter_max = 50, verbose = T)
max(result$match_id, na.rm = T)
resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)

#Changing the order does not change the number of matched sets
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "4-3-2-1",
                    exactMatch = ~var1+var2,
                    iterate = T, niter_max = 50, verbose = T)
max(result$match_id, na.rm = T)
resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)

#If exact match does not exists, not all of the subjects are matched
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean",
                    start = "1-2-3-4",
                    exactMatch = ~var1,
                    iterate = T, niter_max = 50, verbose = T)
max(result$match_id, na.rm = T)
resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)

#Exact match (3) - Starting matched sets is/isn't exactly matched
#-----------------------------------------------------------------
set.seed(123456)
dat <- generateData(c(100,100,100,100),
                    par = list(means = c(0,0,0,.5),
                               sds = c(1,1,3,1)))
dat$var1 <- factor(apply(rmultinom(nrow(dat), size = 1 , prob = c(1/10,2/10,3/10,4/10)), FUN = function(x){which(x==1)},2))
dat$var2 <- factor(rbinom(nrow(dat), size = 1 , prob = c(1/10,9/10)), levels = c(0,1), labels = c("A","B"))

#Initial matched set is not exactly matched
dat$match_id <- NA
for(group in unique(dat$group)) {
  dat$match_id[dat$group %in% group][sample(sum(dat$group %in% group),100)] <- 1:100
}

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    start = dat$match_id,
                    exactMatch = ~var1+var2)

#generate an initial matched set exactly matched and use it as starting point
result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    start = "1-2-3-4",
                    exactMatch = ~var1+var2, iterate = F)

result2 <- polymatch(formulaMatch = group ~ variable, data = dat,
                    start = result$match_id,
                    exactMatch = ~var1+var2, iterate = T)

table(dat$group,dat$var1, dat$var2, dnn = c("group", "var1","var2"))

resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result$match_id)
resultPlot <- plotBalance(resultBalance)

resultBalance <- balance(group ~ variable + var1 + var2,
                         data = dat, match_id = result2$match_id)
resultBalance
resultPlot <- plotBalance(resultBalance)

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


# #Plot distribution of matching variable
# library(ggplot2)
# ggplot(data = dat) + geom_density(aes(x=variable,
#                                       colour = group,
#                                       group = group))
#
# ggplot(data = dat[!is.na(dat$match_id),]) +
#   geom_density(aes(x=variable,
#                    colour = group,
#                    group = group))

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
resultPlot <- plotBalance(resultBalance,ratioVariances = T)

#Check balance (3) - our matched dataset
#-----------------------------------------------
load("C:/Users/giova/Google Drive/PhD/Classes/Individual study - Propensity score/NCH data/bestResult_3wayconstr_allData_1-2.Rdata")
dat <- result3wayConstr$matchedData

dat$Age <- dat$AGE
dat$ISS <- dat$iss
dat$Income.Q1 <- factor(dat$income_1)
dat$Income.Q2 <- factor(dat$income_2)
dat$Income.Q3 <- factor(dat$income_3)
dat$Income.Q4 <- factor(dat$income_4)
dat$Payer.Medicare <- factor(dat$pay1_1)
dat$Payer.Medicaid <- factor(dat$pay1_2)
dat$Payer.PrivInsur <- factor(dat$pay1_3)
dat$Payer.SelfPay <- factor(dat$pay1_4)
dat$Payer.NoCharge <- factor(dat$pay1_5)
dat$Payer.Other <- factor(dat$pay1_6)
dat$Loc.LargeMetroCentral <- factor(dat$nchs_1)
dat$Loc.LargeMetroFringe <- factor(dat$nchs_2)
dat$Loc.MediumMetro <- factor(dat$nchs_3)
dat$Loc.SmallMetro <- factor(dat$nchs_4)
dat$Loc.Micro <- factor(dat$nchs_5)
dat$Loc.Other <- factor(dat$nchs_6)
dat$Chron.Condition <- factor(dat$chronic)
dat$Mult.Injury <- factor(dat$multiple_injury)
dat$Female <- factor(dat$FEMALE)

#E.g. 1 - plot on two columns
resultBalance <- balance(HOSP_TRAUMA ~ Age + Female + ISS + Chron.Condition + Mult.Injury +
                           Income.Q1 + Income.Q2 + Income.Q3 + Income.Q4 +
                           Payer.Medicare + Payer.Medicaid + Payer.PrivInsur + Payer.SelfPay + Payer.NoCharge + Payer.Other +
                           Loc.LargeMetroCentral + Loc.LargeMetroFringe + Loc.MediumMetro + Loc.SmallMetro + Loc.Micro + Loc.Other,
                         data = dat, match_id = dat$indexMatch)
resultBalance

resultPlot <- plotBalance(resultBalance, boxplots = FALSE)
# pdf("C:/Users/natt03/Desktop/triplet matching/analysis/balancePlot.pdf", height = 15, width = 6)
# print(resultPlot[[1]])
# dev.off()

pdf("C:/Users/giova/Google Drive/PhD/Classes/Individual study - Propensity score/Presentation JSM/balancePlot.pdf", height = 10, width = 17.8)
resultPlot$plotStdzDiff + ggplot2::facet_wrap(~variable, dir = "v", strip.position = "left", ncol = 3) +
  ggplot2::labs(title=element_blank()) +
  ggplot2::theme(axis.text=element_text(size=20),
                 axis.title=element_text(size=22))
dev.off()

#E.g. 2 - split balance in two
resultBalance1 <- balance(HOSP_TRAUMA ~ Age + Female + ISS + Chron.Condition + Mult.Injury +
                            Payer.Medicare + Payer.Medicaid + Payer.PrivInsur + Payer.SelfPay + Payer.NoCharge + Payer.Other,
                         data = dat, match_id = dat$indexMatch)
resultBalance2 <- balance(HOSP_TRAUMA ~ Income.Q1 + Income.Q2 + Income.Q3 + Income.Q4 +
                            Loc.LargeMetroCentral + Loc.LargeMetroFringe + Loc.MediumMetro + Loc.SmallMetro + Loc.Micro + Loc.Other,
                         data = dat, match_id = dat$indexMatch)
resultPlot1 <- plotBalance(resultBalance1, boxplots = FALSE)
resultPlot2 <- plotBalance(resultBalance2, boxplots = FALSE)

library(gridExtra)
library(ggplot2)
grid.arrange(resultPlot1[[1]] + labs(title=""),resultPlot2[[1]] + labs(title=""), ncol = 2)


#E.g. 3 - plot with ratio variances
resultBalance <- balance(HOSP_TRAUMA ~ AGE + iss +
                           income_1 + income_2 + income_3 + income_4 +
                           pay1_1 + pay1_2 + pay1_3 + pay1_4 + pay1_5 + pay1_6 +
                           nchs_1 + nchs_2 + nchs_3 + nchs_4 + nchs_5 + nchs_6 +
                           chronic + multiple_injury +
                           FEMALE,
                         data = dat, match_id = dat$indexMatch)
resultPlot <- plotBalance(resultBalance, ratioVariances = TRUE)


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

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean", start = "A-B-C",
                    iterate = T, niter_max = 50, verbose = T)

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean", start = "A-C-B",
                    iterate = T, niter_max = 50, verbose = T)

result <- polymatch(formulaMatch = group ~ variable, data = dat,
                    distance = "euclidean", start = "B-C-A",
                    iterate = T, niter_max = 50, verbose = T)

dat$match_id <- result$match_id

#Optimal solution identified (5-8-8 and 10-11-13)
dat

########################
# For examples in help #
########################

#Generate a datasets with group indicator and four variables:
#- var1, continuous, sampled from normal distributions;
#- var2, continuous, sampled from beta distributions;
#- var3, categorical with 4 levels;
#- var4, binary.
set.seed(12345)
dat <- data.frame(group= c(rep("A",100),rep("B",500),rep("C",500)),
                  var1=c(rnorm(100,mean=0,sd=1),rnorm(500,mean=1,sd=2),rnorm(500,mean=-1,sd=2)),
                  var2=c(rbeta(100,shape1=1,shape2=1),rbeta(500,shape1=2,shape2=1),rbeta(500,shape1=1,shape2=2)),
                  var3=factor(c(rbinom(100,size=3,prob=.4),rbinom(500,size=3,prob=.5),rbinom(500,size=3,prob=.3))),
                  var4=factor(c(rbinom(100,size=1,prob=.5),rbinom(500,size=1,prob=.3),rbinom(500,size=1,prob=.7))))

#Match on propensity score
#-------------------------

#With multiple groups, need a multinomial model for the PS
library(VGAM)
psModel <- vglm(group ~ var1 + var2 + var3 + var4,
                family=multinomial, data=dat)
#Estimated probabilities - 3 for each unit: P(group=A), P(group=B), P(group=C)
probsPS <- predict(psModel, type = "response")
dat$probA <- probsPS[,"A"]
dat$probB <- probsPS[,"B"]
dat$probC <- probsPS[,"C"]
#Estimated logits - 2 for each unit: log(P(group=A)/P(group=C)), log(P(group=B)/P(group=C))
logitPS <- predict(psModel, type = "link")
dat$logit_AvsC <- logitPS[,1]
dat$logit_BvsC <- logitPS[,2]

#Match on logits of PS
resultPs <- polymatch(group ~ logit_AvsC + logit_BvsC, data = dat,
                      distance = "euclidean")
dat$match_id_ps <- resultPs$match_id

#Compare the distributions of propensity score before and after matching
library(ggplot2)
library(tidyr)
library(gridExtra)

#Distribution of propensity score BEFORE matching
distrPsBefore <- ggplot(dat %>%
                          gather(key = probGroup,
                                 value = prob, probA, probB, probC)) +
  geom_density(aes(prob,stat(count),colour=group)) +
  facet_wrap(~factor(probGroup))
#Distribution of propensity score AFTER matching
distrPsAfter <- ggplot(dat %>%
                         drop_na(match_id_ps) %>%
                         gather(key = probGroup,
                                value = prob, probA, probB, probC)) +
  geom_density(aes(prob,stat(count),colour=group)) +
  facet_wrap(~factor(probGroup))
#Single plot with the two distributions
grid.arrange(distrPsBefore +
               labs(title="Distribution of PS before matching"),
             distrPsAfter +
               labs(title="Distribution of PS after matching"),
             nrow = 2)

#Evaluate balance in covariates
tabBalancePs <- balance(group ~ var1 + var2 + var3 + var4,
                        match_id = dat$match_id_ps, data = dat)
tabBalancePs
plotPs <- plotBalance(tabBalancePs, ratioVariances = T)

#Match on covariates
#--------------------

#Match on continuous covariates with exact match on categorical/binary variables
resultCov <- polymatch(group ~ var1 + var2, data = dat,
                       distance = "mahalanobis",
                       exactMatch = ~var3+var4)
dat$match_id_cov <- resultCov$match_id

#Evaluate balance
tabBalanceCov <- balance(group ~ var1 + var2 + var3 + var4,
                         match_id = dat$match_id_cov, data = dat)
tabBalanceCov
plotCov <- plotBalance(tabBalanceCov, ratioVariances = T)

#Compare balance between the two matched samples
#-----------------------------------------------
library(gridExtra)

grid.arrange(plotPs[[1]] +
               labs(title="Stand. Differences - Matching on PS"),
             plotCov[[1]] +
               labs(title="Stand. Differences - Matching on Covariates"),
             ncol = 2)
