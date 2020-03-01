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
#' N=30;
#' L=10;
#' p.in.group =8;
#' p=L * p.in.group;
#' sigma <- sqrt(1);
#' beta.coef <- matrix(0,nrow=2*L,ncol=(p/L)/2)
#' beta.coef[1,] <- c(6,6.4,6.6,8)/2
#' beta.coef[2,] <- c(6,6.4,6.6,8)/2
#' beta.coef[3,] <- c(6,6.6,6.6,8)/2
#' beta.coef[5,] <- c(12.5,12.5,0,0)/2
#' beta.coef <- beta.coef *2
#' p.group <- rep(p/L,L)
#' index.subgroup <- matrix(NA,nrow=L,ncol=p)
#' tmp <- 0
#' for(k in 1:L){
#' if(k==1){
#' index.subgroup[k,1:p.group[k]] <- c(rep(1,(p/L)/2),rep(2,(p/L)/2))
#' } else {
#' ind <- 1:p.group[k] + sum(p.group[(k-1):1])
#' index.subgroup[k,ind] <- c(rep(k+tmp,(p/L)/2),rep(k+tmp+1,(p/L)/2))
#' }
#' tmp <- tmp + 1
#' }
#' out <- data.group(N,p.group,beta.coef,sigma)
#' y <- out$y
#' x <- out$X
#' out_lasso <- sgsl(x,y,type="lasso",index.subgroup = index.subgroup)
#' out_group <- sgsl(x,y,type="group",index.subgroup = index.subgroup,tau=0.94)
#' out_ggroup <- sgsl(x,y,type="ggroup",index.subgroup = index.subgroup,tau=0.94)
#' out_ggroupind <- sgsl(x,y,type="ggroupind",index.subgroup = index.subgroup,tau=0.94)
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
  } else if (type == "group"){
    ## Apply pure group Lasso ##
    out <- pure.grplasso.computations(XX=x,response=y,
                                      index=index,p.group=p.group,tau=tau,
                                      delta.group=delta.group,
                                      standardize=standardize)
  } else if (type == "ggroup"){
    ## Apply the group lasso among subgroups ##
    out <- group.group.lasso.computations(XX=x,response=y,index,index.subgroup,p.group,tau=tau,
                                          delta.group=delta.group,
                                          delta.subgroup=delta.subgroup,
                                          standardize=standardize)

  } else if (type == "ggroupind"){
    ## Apply the lasso with individual features ##
    out <- group.group.indlasso.lasso.computations(XX=x,response=y,index,index.subgroup,p.group,tau=tau,
                                                   delta.group=delta.group,delta.subgroup=delta.subgroup,delta.ind=delta.ind,
                                                   standardize=standardize)

  } else {
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
    lasso.out <- lasso(yy=as.matrix(data.yy[,i]),XX,delta,standardize=standardize)
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




##############################
## PURE GROUP LASSO Approach #
##############################

pure.grplasso.computations <- function(XX,response,index,p.group,tau,
                                       delta.group=2,
                                       standardize=TRUE){

  interest <- matrix(0,nrow=ncol(XX),ncol=ncol(response))
  interest <- as.data.frame(interest)
  rownames(interest) <- colnames(XX)
  colnames(interest) <- colnames(response)

  data.yy <- response


  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- group.lasso(data.yy[,i],XX,index,tau,
                             delta=delta.group,standardize=standardize)

    interest[,i] <- lasso.out
  }
  return(interest)
}




#########################
## GROUP LASSO Approach #
#########################

#' @import stats
#' @import grplasso
group.lasso <- function(yy,XX,index,tau,
                        delta=2,standardize=TRUE){
  y <- yy
  X <- XX

  N <- length(y)

  if(standardize==TRUE){
    y1 <- make.std(y)
    X1 <- apply(X,2,make.std)
  } else {
    y1 <- y
    X1 <- X
  }

  ## Use a multiplicative grid for penalty parameter lambda, starting at maximal lambda value
  lambda <- grplasso::lambdamax(X1,y=y1, index=index,penscale=sqrt,model=LinReg(),
                                center=FALSE,
                                standardize=standardize)  * c(tau^(0:100),0)

  ## Run Group Lasso
  Lasso.out <-  grplasso::grplasso(X1, y = y1, index = index, lambda = lambda,
                                   model = LinReg(),
                                   penscale = sqrt,
                                   control = grpl.control(update.hess = "lambda",
                                                          trace = 0),center=FALSE,
                                   standardize=standardize)

  ## use Cp-like criterion to find best descriptive model
  p = dim(X1)[2]
  s = length(Lasso.out$lambda)
  p.pos = NULL

  RSS = NULL
  for (i in 1:s){
    RSS[i] = sum((y1-Lasso.out$fitted[,i])**2)
    p.pre = Lasso.out$coefficients[,i]
    p.pos = c(p.pos,length(p.pre[abs(p.pre)>1e-10]))
  }

  ## Estimated MSE
  MSE <- stats::sd(as.vector(y1)) * sqrt( N / (N-1) )
  MSE <- MSE^2


  p.min = which.min(RSS/MSE+delta*p.pos)

  ## final best descriptive model
  predict.out <- Lasso.out$coefficients[,p.min]

  ind <- which(abs(predict.out)>1e-10)
  sig.variables <- rep(0,ncol(XX))
  sig.variables[ind] <- 1


  return(sig.variables)
}




######################################################
## PURE GROUP LASSO, GROUP LASSO, and LASSO Approach #
######################################################
## delta.group : delta applied to C_p criterion for group lasso
## delta.subgroup : delta applied to C_p critierian for group lasso among subgroups
## delta.ind    : delta applied to C_p criterion for lasso with individual features

group.group.indlasso.lasso <- function(response,XX,index,index.subgroup,p.group,tau,
                                       delta.group=2,delta.subgroup=2,delta.ind=2,
                                       standardize=TRUE){

  ## group lasso to groups
  group.lasso.out <- group.lasso(response,XX,index,tau,
                                 delta=delta.group,standardize=standardize)

  sig.variables <- group.lasso.out

  main.ind <- which(sig.variables==1)

  ################################################
  ## setting up to apply group Lasso to subgroups #
  #################################################

  tmp.subgroup <- as.vector(t(index.subgroup))
  tmp.subgroup <- tmp.subgroup[!is.na(tmp.subgroup)]
  new.index.subgroup <- tmp.subgroup[main.ind]

  if(sum(main.ind)!=0){
    X.tmp <- XX[,main.ind]
    ##if(length(main.ind)>1){
    ##  X.tmp <- microbes[main.ind,]
    ##} else {
    ##   X.tmp <- t(data.frame(microbes[main.ind,]))
    ##}
    subgroup.lasso.out <- group.lasso(response,X.tmp,new.index.subgroup,tau,
                                      delta=delta.subgroup,standardize=standardize)
    sig.variables[main.ind] <- subgroup.lasso.out
  }

  ######################################################
  ## setting up to apply lasso among variables selected #
  #######################################################

  new.main.ind <- which(sig.variables==1)

  if(sum(new.main.ind)!=0){
    X.tmp <- XX[,new.main.ind]

    ##if(length(new.main.ind)>1){
    ##  X.tmp <- microbes[new.main.ind,]
    ##} else {
    ##   X.tmp <- t(data.frame(microbes[new.main.ind,]))
    ##}
    lasso.out <-  lasso(response,X.tmp,
                        delta=delta.ind)
    sig.variables[new.main.ind] <- lasso.out
  }

  return(sig.variables)

}

group.group.indlasso.lasso.computations <- function(XX,response,index,index.subgroup,p.group,tau,
                                                    delta.group=2,delta.subgroup=2,delta.ind=2,
                                                    standardize=TRUE){
  interest <- matrix(0,nrow=ncol(XX),ncol=ncol(response))
  interest <- as.data.frame(interest)
  rownames(interest) <- colnames(XX)
  colnames(interest) <- colnames(response)

  data.yy <- response

  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- group.group.indlasso.lasso(data.yy[,i],XX,index,index.subgroup,p.group,tau,
                                            delta.group,delta.subgroup,delta.ind,
                                            standardize=standardize)

    interest[,i] <- lasso.out
  }
  return(interest)
}



############################################
## PURE GROUP LASSO & GROUP LASSO Approach #
############################################

## delta.group : delta applied to C_p criterion for group lasso
## delta.subgroup : delta applied to C_p critierian for group lasso among subgroups

group.group.lasso <- function(yy,XX,index,index.subgroup,p.group,tau,
                              delta.group=2,delta.subgroup=2,
                              standardize=TRUE){

  ## group lasso to groups
  group.lasso.out <- group.lasso(yy,XX,index,tau,
                                 delta=delta.group,standardize=standardize)

  sig.variables <- group.lasso.out

  main.ind <- which(sig.variables==1)

  ## setting up to apply group Lasso to subgroups
  tmp.subgroup <- as.vector(t(index.subgroup))
  tmp.subgroup <- tmp.subgroup[!is.na(tmp.subgroup)]
  new.index.subgroup <- tmp.subgroup[main.ind]

  if(sum(main.ind)!=0){
    X.tmp <- XX[,main.ind]

    subgroup.lasso.out <- group.lasso(yy,X.tmp,new.index.subgroup,tau,
                                      delta=delta.subgroup,standardize=standardize)
    sig.variables[main.ind] <- subgroup.lasso.out
  }

  return(sig.variables)

}

group.group.lasso.computations <- function(XX,response,index,index.subgroup,p.group,tau,
                                           delta.group=2,
                                           delta.subgroup=2,standardize=TRUE){

  interest <- matrix(0,nrow=ncol(XX),ncol=ncol(response))
  interest <- as.data.frame(interest)
  rownames(interest) <- colnames(XX)
  colnames(interest) <- colnames(response)

  data.yy <- response

  for(i in 1:ncol(interest)){
    ##print(i)
    lasso.out <- group.group.lasso(data.yy[,i],XX,
                                   index,index.subgroup,p.group,tau,
                                   delta.group,delta.subgroup,
                                   standardize=standardize)

    interest[,i] <- lasso.out
  }
  return(interest)
}



