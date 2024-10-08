# Polymatching: Matching in Designs with Multiple Treatment Groups

## Description

The package implements the conditionally optimal matching algorithm, which can be used to generate matched samples in designs with multiple treatment groups.

## Details

Currently, the algorithm can be applied to datasets with 3 to 6 groups and generates matched samples with one subject per group. The package provides functions to generate the matched sample and to evaluate the balance in key covariates.

## Generating the Matched Sample

The function implementing the matching algorithm is `polymatch`. The algorithm is iterative and needs a matched sample with one subject per group as starting point. This matched sample can be automatically generated by `polymatch` or can be provided by the user. The algorithm iteratively explores possible reductions in the total distance of the matched sample.

## Evaluating Balance in Covariates

Balance in key covariates can be evaluated with the function `balance`. Given a matched sample and a set of covariates of interest, the function computes the standardized differences and the ratio of the variances for each pair of treatment groups in the study design. For 3, 4, 5 and 6 groups, there are 3, 6, 10 and 15 pairs of groups and the balance is evaluated before and after matching. The result of balance can be graphically represented with `plotBalance`.


## Installation

You can install the package with the function `install_github` of the package `devtools`.

```
library(devtools)
install_github("gnattino/polymatching")
```
