#' Fit the sparse group-subgroup lasso (SGSL)
#'
#' @param x p by N matrix of predictors (N: sample size, p: number of predictors)
#' @param y 1 by N matrix of response variable
#' @param type One of "lasso" (for the standard lasso), "group" (for the group lasso), "ggroup" (for the group lasso among subgroups), "ggroupind" (for the lasso with individual features), "sgsl" (for the sparse-group-subgroup lasso) or "groupsgl (for the sparse group lasso at subgroup level).
#' @param alpha1 REgularization parameter.
#' @param alpha2 regularization parameter.
#' @param alpha3 regularization parameter.
#' @param cv.criterion logical indicator.
#' @param nlam number of lambda values.
#' @param index.subgroup index for subgroups
#' @param tau multiplier for using a multiplicative grid for penalty parameter lambda, starting at maximal lambda value
#' @param delta Among the lasso solution path, the best descriptive model is the one which minimizes the loss function: (residual sum of squares)/(estimator of the model error variance) - (sample size) + delta*(number of predictors in the selected model). If delta = 2, this loss function is Mallows' Cp.
#' @param delta.group delta applied to C_p criterion for group lasso
#' @param delta.subgroup delta applied to C_p critierian for group lasso among subgroups
#' @param delta.ind delta applied to C_p criterion for lasso with individual features
#' @param standardize logical. TRUE for standardizing the data.
#'
#' @return out: indicators of the selected predictors. 1 for selected predictors and 0 for not selected predictors
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(360), nrow=12)
#' y <- 0.5*x[1,] + 0.5*x[2,] + 1.0*x[4,] + matrix(rnorm(30), nrow=1)
#' index.subgroup <- matrix(NA,nrow=3,ncol=12)
#' index.subgroup[1,1:2]=1; index.subgroup[1,3:4]=2
#' index.subgroup[2,5:6]=3; index.subgroup[2,7:8]=4
#' index.subgroup[3,9:10]=5; index.subgroup[3,11:12]=6
#' out_lasso <- sgsl(x,y,type="lasso",index.subgroup = index.subgroup)
#' ##out_group <- sgsl(x,y,type="group",index.subgroup = index.subgroup,tau=0.94)
#' ##out_ggroup <- sgsl(x,y,type="ggroup",index.subgroup = index.subgroup,tau=0.94)
#' ##out_ggroupind <- sgsl(x,y,type="ggroupind",index.subgroup = index.subgroup,tau=0.94)
#'
sgsl <- function(x,y,type=c("lasso", "group", "ggroup", "ggroupind", "sgsl", "groupsgl")[1],
                 index.subgroup,
                 tau=0.94,
                 alpha1=0.05,
                 alpha2=0.2,
                 alpha3= 0.1,
                 cv.criterion=FALSE,
                 nlam=100,
                 delta=2,delta.group=2,delta.subgroup=2,delta.ind=2,
                 standardize=TRUE){

  tools <- form.tools(index.subgroup)
  index <- tools$index
  p.group <- tools$p.group
  p.subgroup <- tools$p.subgroup
  x <- as.data.frame(x)
  y <- as.data.frame(y)

  ## allocate names to X
  if(is.null(colnames(x))){
    colnames(x) <- paste0("X_",1:ncol(x))
  }

  if(is.null(colnames(y))){
    colnames(y) <- paste0("y",1:ncol(y))
  }

  ## Methods ##
  if (type == "lasso"){
    ## Apply Lasso ##
    out <- lasso.computations(XX=x,response=y,delta=delta,standardize=standardize)
  }  else {
    print("type should be either lasso, group, ggroup, ggroupind, sgsl or groupsgl")
  }

  return(out)

}


###################################
## Functions for organizing data ##
###################################

## function to count number of non-NA elements in a vector

number.non.na <- function(x){
  return(sum(!is.na(x)))
}



## function to form p.group,p.subgroup, and index variables
## p.group : number of covariates in each group
## p.subgroup : matrix of number of covariates in each subgroup (each row corresponds to a group)
## index : group assignments

form.tools <- function(index.subgroup){
  p.group <- apply(index.subgroup,1,number.non.na)
  p.subgroup <- as.numeric(table(index.subgroup))
  index <- rep(1:length(p.group),p.group)
  list(p.group=p.group,p.subgroup=p.subgroup,index=index)
}





####################
## LASSO Approach ##
####################


lasso.computations <- function(XX,response,delta,standardize){

  interest <- matrix(0,nrow=ncol(XX),ncol=ncol(response))
  interest <- as.data.frame(interest)
  rownames(interest) <- colnames(XX)
  colnames(interest) <- colnames(response)

  data.yy <- response

  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- lasso(yy=data.yy[,i],XX,delta,standardize=standardize)
    interest[,i] <- lasso.out
  }
  return(interest)
}




## function to standardize a vector x

#' @import stats
make.std <- function(x){
  N <- length(x)
  ( x-mean(x) ) / ( stats::sd(as.vector(x)) * sqrt( N / (N-1) ) )
}

## function to center a vector x

make.center <- function(x){
  return(x-mean(x))
}

#' @import stats
#' @import lars
lasso <- function(yy,XX,delta=2,standardize=TRUE){
  #####################
  ## Set up the data ##
  #####################

  y <- yy
  X <- XX

  N <- length(y)

  if(standardize==TRUE){
    y1 <- make.std(y)

    ## Standardized design matrix X
    X1 <- apply(X,2,make.std)
  } else {
    y1 <- y
    X1 <- X
  }

  ## adjust use.Gram
  if(ncol(X1)>500){
    use.Gram <- FALSE
  } else {
    use.Gram <- TRUE
  }

  ## Run Lasso
  Lasso.out <-  lars::lars(X1, y1, type = c("lasso"),
                           trace = FALSE, normalize = FALSE, intercept = FALSE,
                           use.Gram=use.Gram)

  ## use Cp-like criterion to find best descriptive model
  p = dim(X1)[2]
  s = length(Lasso.out$df)
  p.pos = NULL

  RSS = NULL
  for (i in 1:s){
    RSS[i] = sum((y1-lars::predict.lars(Lasso.out, X1, s=i, type = c("fit"))$fit)**2)
    p.pre = lars::predict.lars(Lasso.out, X1, s=i, type = c("coefficients"))$coefficients
    p.pos = c(p.pos,length(p.pre[abs(p.pre)>1e-10]))
  }

  # Estimated MSE
  MSE <- stats::sd(as.vector(y1)) * sqrt( N / (N-1) )
  MSE <- MSE^2

  p.min = which.min(RSS/MSE+delta*p.pos)

  ## final best descriptive model
  predict.out <- lars::predict.lars(Lasso.out, X1, s=p.min, type = c("coefficients"))

  ind <- which(abs(predict.out$coefficients)>1e-10)
  sig.variables <- rep(0,ncol(XX))
  sig.variables[ind] <- 1

  return(sig.variables)
}


