% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polymatch.R
\name{polymatch}
\alias{polymatch}
\title{Polymatching}
\usage{
polymatch(
  formulaMatch,
  start = "small.to.large",
  data,
  distance = "euclidean",
  exactMatch = NULL,
  vectorK = NULL,
  iterate = TRUE,
  niter_max = 50,
  withinGroupDist = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{formulaMatch}{Formula with form \code{group ~ x_1 + ... + x_p}, where \code{group} is the name of the variable
identifying the treatment groups/exposures and \code{x_1},...,\code{x_p} are the matching variables.}

\item{start}{An object specifying the starting point of the iterative algorithm. Three types of input are accepted:
\itemize{
  \item \code{start="small.to.large"} (default): the starting matched set is generated by matching groups
from the smallest to the largest.
  \item Users can specify the order to be used to match groups for the starting sample.
For example, if there are four groups with labels "A","B","C" and "D", \code{start="D-B-A-C"} generates the starting sample
by matching groups "D" and "B", then units from "A" to the "D"-"B"pairs, then units from "C" to the "D"-"B"-"A" triplets.
  \item Users can provide the starting matched set and the algorithm will explore possible reductions in the total
distance. In this case, \code{start} must be a vector with the IDs of the matched sets, i.e.,
a vector with length equal to the number of rows of \code{data} where
matched subjects are flagged with the same value and non-matched subjects have value \code{NA}.
}}

\item{data}{The \code{data.frame} object with the data.}

\item{distance}{String specifying whether the distance between pairs of observations should be computed with the
Euclidean (\code{"euclidean"}, default) or Mahalanobis (\code{"mahalanobis"}) distance. See section 'Details' for further information.}

\item{exactMatch}{Formula with form \code{~ z_1 + ... + z_k}, where \code{z_1},...,\code{z_k} must
be factor variables. Subjects are exactly matched on \code{z_1},...,\code{z_k}, i.e., matched
within levels of these variables.}

\item{vectorK}{A named vector with the number of subjects from each group in each matched set. The names of the vector must be 
the labels of the groups, i.e., the levels of the variable identifying the treatment groups/exposures. 
For example, in case of four groups with labels "A","B","C" and "D" and assuming that the desired design is 1:2:3:3 
(1 subject from A, 2 from B, 3 from C and 3 from D in each matched set), the parameter should be set to
\code{vectorK =  c("A" = 1, "B" = 2, "C" = 3, "D" = 3)}. By default, the generated matched design includes 1 subject per group in each
matched set, i.e, a 1:1: ... :1 matched design.}

\item{iterate}{Boolean specifying whether iterations should be done (\code{iterate=TRUE}, default) or not (\code{iterate=FALSE}).}

\item{niter_max}{Maximum number of iterations. Default is 50.}

\item{withinGroupDist}{Boolean specifying whether the distances within the same treatment/exposure group should be considered in the 
total distance. For example, in a 1:2:3 matched design among the groups A, B and C, the parameters controls whether the distance 
between the two subjects in B and the three pairwise distances among the subjects in C should be counted in the total distance. 
The default value is \code{TRUE}.}

\item{verbose}{Boolean: should text be printed in the console? Default is \code{TRUE}.}
}
\value{
A list containing the following components:
\describe{
  \item{match_id}{A numeric vector identifying the matched sets---matched units have the same identifier.}
  \item{total_distance}{Total distance of the returned matched sample.}
  \item{total_distance_start}{Total distance at the starting point.}
}
}
\description{
\code{polymatch} generates matched samples in designs with up to 10 groups.
}
\details{
The function implements the conditionally optimal matching algorithm, which iteratively uses
two-group optimal matching steps to generate matched samples with small total distance. In the current implementation,
it is possible to generate matched samples with multiple subjects per group, with the matching ratio being 
specified by the \code{vectorK} parameter.

The steps of the algorithm are described with the following example. Consider a 4-group design with
groups labels "A", "B", "C" and "D" and a 1:1:1:1 matching ratio. The algorithm requires a set of quadruplets as starting point. 
The argument \code{start} defines the approach to be used to
generate such a starting point. \code{polymatch} generates the starting point by sequentially using optimal two-group matching.
In the default setting (\code{start="small.to.large"}), the steps are:
\enumerate{
 \item optimally match the two smallest groups;
 \item optimally match the third smallest group to the pairs generated in the first step;
 \item optimally match the last group to the triplets generated in the second step.
}
Notably, we can use the optimal two-group algorithm in steps 2) and 3) because they are
two-dimensional problems: the elements of one group on one hand, fixed matched sets on the other hand. The order of the
groups to be considered when generating the starting point can be user-specified (e.g., \code{start="D-B-A-C"}).
In alternative, the user can provide a matched set that will be used as starting point.

Given the starting matched set, the algorithm iteratively explores possible reductions in the total distance (if \code{iterate="TRUE"}),
by sequentially relaxing the connection to each group and rematching units of that group. In our example:
\enumerate{
 \item rematch "B-C-D" triplets within the starting quadruplets to units in group "A";
 \item rematch "A-C-D" triplets within the starting quadruplets to units in group "B";
 \item rematch "A-B-D" triplets within the starting quadruplets to units in group "C";
 \item rematch "A-B-C" triplets within the starting quadruplets to units in group "D".
}
If none of the sets of quadruplets generated in 1)-4) has smaller total distance than the starting point, the algorihm stops.
Otherwise, the set of quadruplets with smallest distance is seleceted and the process iterated, until no reduction in the total
distance is found or the number of maximum iterations is reached (\code{niter_max=50} by default).

The total distance is defined as the sum of all the within-matched-set distances. The within-matched-set distance is defined as the
sum of the pairwise distances between pairs of units in the matched set. The type of distance is specified with the \code{distance}
argument. The current implementation supports Euclidean (\code{distance="euclidean"}) and Mahalanobis (\code{distance="mahalanobis"})
distances. In particular, for the Mahalanobis distance, the covariance matrix is defined only once on the full dataset.
}
\examples{
#Generate a datasets with group indicator and four variables:
#- var1, continuous, sampled from normal distributions;
#- var2, continuous, sampled from beta distributions;
#- var3, categorical with 4 levels;
#- var4, binary.
set.seed(1234567)
dat <- data.frame(group = c(rep("A",20),rep("B",60),rep("C",60)),
                  var1 = c(rnorm(20,mean=0,sd=1),
                           rnorm(60,mean=1,sd=2),
                           rnorm(60,mean=-1,sd=2)),
                  var2 = c(rbeta(20,shape1=1,shape2=1),
                           rbeta(60,shape1=2,shape2=1),
                           rbeta(60,shape1=1,shape2=2)),
                  var3 = factor(c(rbinom(20,size=3,prob=.4),
                                  rbinom(60,size=3,prob=.5),
                                  rbinom(60,size=3,prob=.3))),
                  var4 = factor(c(rbinom(20,size=1,prob=.5),
                                  rbinom(60,size=1,prob=.3),
                                  rbinom(60,size=1,prob=.7))))

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
resultPs1 <- polymatch(group ~ logit_AvsC + logit_BvsC, data = dat,
                    distance = "euclidean")
dat$match_id_ps1 <- resultPs1$match_id

#Match on probabilities
resultPs2 <- polymatch(group ~ probA + probB + probC, data = dat,
                       distance = "euclidean")
dat$match_id_ps2 <- resultPs2$match_id

#Match on covariates
#--------------------

#Match on continuous covariates only
resultCov1 <- polymatch(group ~ var1 + var2, data = dat,
                        distance = "mahalanobis")
dat$match_id_cov1 <- resultCov1$match_id

#Match on continuous covariates with exact match on categorical/binary variables
resultCov2 <- polymatch(group ~ var1 + var2, data = dat,
                        distance = "mahalanobis",
                        exactMatch = ~var3+var4)
dat$match_id_cov2 <- resultCov2$match_id

}
\seealso{
\code{\link{balance}} and \code{\link{plotBalance}} to summarize the
balance in the covariates.
}
