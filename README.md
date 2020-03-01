---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->


# sgslasso

<!-- badges: start -->
<!-- badges: end -->

This package implements the sparse group-subgroup lasso which  identifies important regressor groups, subgroups and individual variables. 
Results from the sparse group subgroup lasso (SGSL) are compared against combinations of group lasso and lasso. The corresponding references are:

Garcia, T.P., Müller, S., Carroll, R.J., and Walzem, R.L. (2013).
Identification of important regressor groups, subgroups, and individuals via regularization methods: application to gut
microbiome data. Bioinformatics, DOI: 10.1093/bioinformatics/btt608.

## Installation

You can install the released version of sgslasso from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sgslasso")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tpgarcia/sgslasso", force=TRUE)
```
## Example

This is a simple example of using SGSL:

```r
library(sgslasso)

## Generate the data
set.seed(1)
N=30;
L=10;
p.in.group =8;
p=L * p.in.group;
sigma <- sqrt(1);
beta.coef <- matrix(0,nrow=2*L,ncol=(p/L)/2)
beta.coef[1,] <- c(6,6.4,6.6,8)/2
beta.coef[2,] <- c(6,6.4,6.6,8)/2
beta.coef[3,] <- c(6,6.6,6.6,8)/2
beta.coef[5,] <- c(12.5,12.5,0,0)/2
beta.coef <- beta.coef *2
p.group <- rep(p/L,L)
index.subgroup <- matrix(NA,nrow=L,ncol=p)
tmp <- 0
for(k in 1:L){
  if(k==1){
    index.subgroup[k,1:p.group[k]] <- c(rep(1,(p/L)/2),rep(2,(p/L)/2))
  } else {
    ind <- 1:p.group[k] + sum(p.group[(k-1):1])
    index.subgroup[k,ind] <- c(rep(k+tmp,(p/L)/2),rep(k+tmp+1,(p/L)/2))
  }
  tmp <- tmp + 1
}
out <- data.group(N,p.group,beta.coef,sigma)
y <- out$y
x <- out$X

## Lasso variable selection
out_lasso <- sgsl(x,y,type="lasso",index.subgroup = index.subgroup)

## Group lasso variable selection among the groups only
##out_group <- sgsl(x,y,type="group",index.subgroup = index.subgroup,tau=0.94)

## Group lasso variable selection among the groups and subgroups
##out_ggroup <- sgsl(x,y,type="ggroup",index.subgroup = index.subgroup,tau=0.94)

## Group lasso variable selection among the groups and subgroups, and lasso among variables within each subgroup
##out_ggroupind <- sgsl(x,y,type="ggroupind",index.subgroup = index.subgroup,tau=0.94)

## Sparse group subgroup lasso
out_sgsl <- sgsl(x,y,type="sgsl",index.subgroup=index.subgroup,tau=0.94)
```
