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
#' out_sgsl <- sgsl(x,y,type="sgsl",index.subgroup=index.subgroup,tau=0.94)
sgsl <- function(x,y,type=c("lasso", "group", "ggroup", "ggroupind", "sgsl")[1],
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
  } else if (type == "sgsl"){


    out <- sparse.group.subgroup.computations(XX=x,response=y,group.index=index,subgroup.index=index.subgroup,
                                              tau=tau,alpha1,alpha2,alpha3,nlam,
                                              lambdas=NULL,lambda.accuracy=1e-4,
                                              delta.group=delta.group,
                                              cv.criterion=cv.criterion,
                                              nfold=10,alphas.cv.range=seq(0.1,0.95,by=0.05))
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


#' Title
#'
#' @param XX
#' @param response
#' @param delta
#' @param standardize
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param XX
#' @param response
#' @param index
#' @param index.subgroup
#' @param p.group
#' @param tau
#' @param delta.group
#' @param delta.subgroup
#' @param delta.ind
#' @param standardize
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param XX
#' @param response
#' @param index
#' @param index.subgroup
#' @param p.group
#' @param tau
#' @param delta.group
#' @param delta.subgroup
#' @param standardize
#'
#' @return
#' @export
#'
#' @examples
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











################################
## SPARSE GROUP SUBGROUP LASSO #
################################

#' @importFrom SGL predictSGL
#' @import stats
sparse.group.subgroup <- function(yy,XX,group.index,subgroup.index,tau,alpha1=0.45,alpha2=0.45,
                                  alpha3=1-alpha1-alpha2,
                                  nlam=100,lambdas=NULL,lambda.accuracy=1e-4,
                                  delta=2,standardize=TRUE){
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

  ## Put data in a list
  data.list <- list(x=X1,y=y1)


  ## Run sparse group subgroup lasso
  Lasso.out <-  subgroup.SGL(data.list, group.index=group.index,
                             subgroup.index=subgroup.index,
                             type = "linear",nlam=nlam,standardize=FALSE,
                             alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,
                             lambdas=lambdas,
                             lambda.accuracy=lambda.accuracy)
  if(ncol(X1)==1){
    Lasso.out$beta <- t(as.matrix(Lasso.out$beta))
  }


  ## Samuel's corrections 11/9/2011
  ## use Cp-like criterion to find best descriptive model
  p = dim(X1)[2]
  s = length(Lasso.out$lambdas)
  p.pos = NULL

  RSS = NULL
  for (i in 1:s){
    RSS[i] = sum((y1-SGL::predictSGL(Lasso.out, X1, lam=i))**2)
    p.pre = Lasso.out$beta[,i]
    p.pos = c(p.pos,length(p.pre[abs(p.pre)>1e-10]))
  }

  ## Estimated MSE
  MSE <- stats::sd(as.vector(y1)) * sqrt( N / (N-1) )
  MSE <- MSE^2

  p.min = which.min(RSS/MSE+delta*p.pos)
  ##print("subgroup SGL")
  ##print(Lasso.out$lambdas[p.min])

  ## final best descriptive model
  predict.out <- Lasso.out$beta[,p.min]

  ind <- which(abs(predict.out)>1e-10)
  sig.variables <- rep(0,ncol(XX))
  sig.variables[ind] <- 1

  return(list(sig.variables=sig.variables,predict.out=predict.out))
}

## cross-validation code to get one set of alphas for sparse.group.subgroup

cv.sparse.group.subgroup.alphas <- function(yy,XX,group.index,subgroup.index,tau,alphas=NULL,
                                            nfold=10,ratio=1,
                                            nlam=100,lambdas=NULL,lambda.accuracy=1e-4,
                                            delta=2,standardize=TRUE){

  ## Matrix to store results
  residmat <- matrix(NA, nrow(alphas), nfold)

  ## do transpose to get correct dimension
  y <- yy
  X <- XX

  ## Randomly partition the data
  all.folds <- cv.folds(length(y),nfold)

  for(a in 1:nrow(alphas)){
    for(j in seq(nfold)){
      ## data to omit
      omit <- all.folds[[j]]

      yy.fold <- t(y[-omit])
      XX.fold <- t(X[-omit,,drop=FALSE])
      beta.omit <- sparse.group.subgroup(yy.fold,XX.fold,group.index,subgroup.index,tau,
                                         alpha1=as.numeric(alphas[a,1]),
                                         alpha2=as.numeric(alphas[a,2]),
                                         alpha3=as.numeric(alphas[a,3]),
                                         nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                         delta=delta,standardize=standardize)

      ## Find final fit with data omitted
      fit <- X[omit,,drop=FALSE] %*% beta.omit$predict.out

      ## Store residual
      residmat[a,j] <- apply((y[omit]-fit)^2,2,sum)
    }
  }

  cv <- apply(residmat,1,mean)

  ## Check which alpha's lead to min(cv)
  alpha.ind <- which(cv==min(cv,na.rm=TRUE))  ## Ignore NA values
  alpha.opt <- alphas[alpha.ind,]

  ## Is alpha.opt unique? If not, take average.
  if(length(alpha.ind)>1){
    alpha.opt <- apply(alpha.opt,2,mean)
  }

  return(alpha.opt)
}

make.alphas <- function(alphas.cv=NULL,ratio=1){
  ## Index of fixed alpha
  fixed <- which(lapply(alphas.cv,length)==1)

  ## Index of NON-varying alpha
  other <- which(lapply(alphas.cv,is.null)==TRUE)

  ## Index of varying alpha
  vary <- setdiff(1:length(alphas.cv),c(fixed,other))

  alphas <- matrix(0,nrow=length(alphas.cv[[vary]]),ncol=3)
  colnames(alphas) <- c("alpha1","alpha2","alpha3")

  alphas[,vary] <- alphas.cv[[vary]]
  alphas[,other] <- ratio * alphas.cv[[vary]]
  alphas[,fixed] <- alphas.cv[[fixed]]

  ## Only keep rows such that sum(alphas)<=1
  alphas <- good.alphas(alphas)

  return(alphas)
}

## function to keep those alphas such that sum(alphas) <=1

good.alphas <- function(alphas){
  ## Only keep rows such that sum(alphas)<=1
  ind <- which(apply(alphas,1,sum)==1)
  alphas <- alphas[ind,]

  return(alphas)
}

cv.sparse.group.subgroup <- function(yy,XX,group.index,subgroup.index,tau,
                                     alphas.cv.range=c(seq(0.01,0.1,by=0.01),seq(0.15,0.95,by=0.05)),
                                     nfold=10,nlam=100,lambdas=NULL,lambda.accuracy=1e-4,
                                     delta=2,standardize=TRUE,sam.implement=FALSE){

  if(sam.implement==FALSE){
    alphas.cv <- list(alpha1.cv = alphas.cv.range, alpha2.cv = alphas.cv.range, alpha3.cv = alphas.cv.range)
    alphas <- do.call(expand.grid,alphas.cv)
    alphas <- good.alphas(alphas)

    cv.out <- cv.sparse.group.subgroup.alphas(yy,XX,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              standardize=standardize)
    alpha.opt <- as.numeric(cv.out$alpha.opt)

  } else {


    ###########################################################
    ## First application: alpha3=1/3, alpha1=ratio * alpha2, ratio=1 ##
    ###########################################################

    alphas.cv <- list(alpha1.cv = NULL, alpha2.cv = alphas.cv.range, alpha3.cv = 1/3)
    ratio <- 1
    alphas <- make.alphas(alphas.cv,ratio)

    cv.out <- cv.sparse.group.subgroup.alphas(yy,XX,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              standardize=standardize)
    alpha.opt <- cv.out$alpha.opt

    #############################################################################################
    ## Second application: alpha2= alphas[2], alpha1=ratio * alpha3, ratio=alphas[1]/alphas[3] ##
    #############################################################################################

    alphas.cv <- list(alpha1.cv = NULL, alpha2.cv = alpha.opt[2], alpha3.cv = alphas.cv.range)
    ratio <- alpha.opt[1] / alpha.opt[3]
    alphas <- make.alphas(alphas.cv,ratio)

    cv.out <- cv.sparse.group.subgroup.alphas(yy,XX,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              delta=delta,standardize=standardize)
    alpha.opt <- cv.out$alpha.opt


    ############################################################################################
    ## Third application: alpha1= alphas[1], alpha2=ratio * alpha3, ratio=alphas[2]/alphas[3] ##
    ############################################################################################

    alphas.cv <- list(alpha1.cv = alpha.opt[1], alpha2.cv = NULL, alpha3.cv = alphas.cv.range)
    ratio <- alpha.opt[2] / alpha.opt[3]
    alphas <- make.alphas(alphas.cv,ratio)

    cv.out <- cv.sparse.group.subgroup.alphas(yy,XX,group.index,subgroup.index,
                                              tau,alphas=alphas,
                                              nfold=nfold,ratio=ratio,
                                              nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                              delta=delta,standardize=standardize)
    alpha.opt <- cv.out$alpha.opt
  }
  ############################################################
  ## Re-run sparse group-subgroup with optimal alpha values ##
  ############################################################
  Lasso.out <-  sparse.group.subgroup(yy,XX,group.index,subgroup.index,tau,
                                      alpha1=as.numeric(alpha.opt[1]),
                                      alpha2=as.numeric(alpha.opt[2]),
                                      alpha3=as.numeric(alpha.opt[3]),
                                      nlam=nlam,lambdas=lambdas,lambda.accuracy=lambda.accuracy,
                                      delta=delta,standardize=standardize)

  list(sig.variables = Lasso.out$sig.variables, alphas=alpha.opt)
}

#' Title
#'
#' @param XX
#' @param response
#' @param group.index
#' @param subgroup.index
#' @param tau
#' @param alpha1
#' @param alpha2
#' @param alpha3
#' @param nlam
#' @param lambdas
#' @param lambda.accuracy
#' @param delta.group
#' @param format.data
#' @param cv.criterion
#' @param nfold
#' @param alphas.cv.range
#'
#' @return
#' @export
#'
#' @examples
sparse.group.subgroup.computations <- function(XX,response,group.index,subgroup.index,tau,alpha1,alpha2,
                                               alpha3,nlam,lambdas,
                                               lambda.accuracy,
                                               delta.group=2,format.data=TRUE,
                                               cv.criterion=FALSE,nfold=10,alphas.cv.range=seq(0.1,0.95,by=0.05)){

  interest <- matrix(0,nrow=ncol(XX),ncol=ncol(response))
  interest <- as.data.frame(interest)
  rownames(interest) <- colnames(XX)
  colnames(interest) <- colnames(response)

  data.yy <- response


  ## Store alpha values
  alpha.out <- matrix(0,nrow=3,ncol=ncol(response))

  for(i in 1:ncol(interest)){
    ##print(i)
    if(cv.criterion==TRUE){

      lasso.out <- cv.sparse.group.subgroup(data.yy[,i],XX,group.index,subgroup.index,tau,
                                            alphas.cv.range=alphas.cv.range,
                                            nfold=nfold,nlam=nlam,lambdas=lambdas,
                                            lambda.accuracy=lambda.accuracy,
                                            delta=delta.group)
      alpha.out[,i] <- lasso.out

    } else {

      lasso.out <- sparse.group.subgroup(data.yy[,i],XX,group.index,subgroup.index,
                                         tau,alpha1,alpha2,alpha3,nlam=nlam,lambdas=lambdas,
                                         lambda.accuracy=lambda.accuracy,
                                         delta=delta.group)
    }
    interest[,i] <- lasso.out$sig.variables
  }

  return(list(interest=interest,alpha.out=alpha.out))
}








#######################################################
# R function to implement Sparse Group-Subgroup Lasso #
#######################################################

#############
# Functions #
#############

################
# subgroup.SGL #
################

# data : list that contains design matrix (x) and response vector (y)
# group.index: vector indicating which covariates belong to which group
# subgroup.index: matrix indicating which covariates belong to which subgroup within a group,
#		Each row of the subgroup.index corresponds to a group, and entries determine which subgroup the element belongs where
#		some an entry of "NA" means that there is no subgroup association. We ASSUME row 1 of subgroup.index corresponds to elements
#               in group 1, row 2 for elements in group 2, etc.
# type : linear regression model
# maxit: maximum number of iterations in algorithm
# thresh: threshold to determine if beta converges
# min.frac: used in determining appropriate lambda values; this is the minimum value of penalty parameter (lambda_max) as a  fraction of maximum value
# nlam : number of lambda values used
# gamma : fitting parameter used in backtracking (inner loop of subgroup)
# standardize: indicator to standardize covariates
# verbose : indicator to display results as algorithm progresses
# step : fitting parameter used for initial backtracking step size
# reset : fitting parameter used in nesterov momentum; it is the number of iterations before momentum is reset
# alpha1 : alpha_1 value
# alpha2 : alpha_2 value
# lambdas : lambda values used; if NULL, program selects best set of lambda values
# lambda.accuracy : accuracy of lambda.max calculation

subgroup.SGL <- function(data, group.index, subgroup.index,
                         type = "linear", maxit = 1000,
                         thresh = 0.001, min.frac = 0.1, nlam = 20,
                         gamma = 0.8, standardize = TRUE, verbose = FALSE, step = 1, reset = 10,
                         alpha1 = 0.45, alpha2=0.45, alpha3=1-alpha1-alpha2, lambdas = NULL, lambda.accuracy=1e-4){

  X.transform <- NULL

  if(standardize == TRUE){
    X <- data$x
    means <- apply(X,2,mean)
    X <- t(t(X) - means)
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    X <- t(t(X) / var)
    data$x <- X
    X.transform <- list(X.scale = var, X.means = means)
  }

  if(type == "linear"){
    if(standardize == TRUE){
      intercept <- mean(data$y)
      data$y <- data$y - intercept
    }
    Sol <- oneDimNew(data, group.index,subgroup.index, thresh,  nlam = nlam, lambdas = lambdas,
                     inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh,
                     gamma = gamma, step = step, reset = reset, alpha1 = alpha1, alpha2=alpha2, alpha3=alpha3,
                     min.frac = min.frac, verbose = verbose, lambda.accuracy = lambda.accuracy)

    if(standardize == TRUE){
      Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, intercept = intercept, X.transform = X.transform)
    }
    if(standardize == FALSE){
      Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, X.transform = X.transform)
    }
  }

  class(Sol) = "subgroup.SGL"
  return(Sol)
}





##############
# oneDimNew  #
##############

# data : list that contains design matrix (x) and response vector (y)
# group.index: vector indicating which covariates belong to which group
# subgroup.index: matrix indicating which covariates belong to which subgroup within a group,
#		Each row of the subgroup.index corresponds to a group, and entries determine which subgroup the element belongs where
#		some an entry of "NA" means that there is no subgroup association.
# thresh: threshold to determine if beta converges
# nlam : number of lambda values used
# lambdas : lambda values used; if NULL, program selects best set of lambda values
# beta.naught : initial beta values, these must be zero!!!
# inner.iter : maximum number of iterations for inner loop (set to maxit)
# outer.iter : maximum number of iterations for outer loop (set to maxit)
# outer.thresh : threshold to determine convergence in iterations for outer loop
# gamma : fitting parameter used in backtracking (inner loop of subgroup)
# step : fitting parameter used for initial backtracking step size
# reset : fitting parameter used in nesterov momentum; it is the number of iterations before momentum is reset
# alpha1 : alpha_1 value
# alpha2 : alpha_2 value
# min.frac: used in determining appropriate lambda values; this is the minimum value of penalty parameter (lambda_max) as a  fraction of maximum value
# verbose : indicator to display results as algorithm progresses
oneDimNew <-
  function(data, group.index, subgroup.index, thresh = 0.0001, nlam = 20, lambdas = NULL,
           beta.naught = rep(0,ncol(data$x)), inner.iter = 100, outer.iter = 100, outer.thresh = 0.0001,
           gamma = 0.8, step = 1, reset = 10, alpha1 = 0.45, alpha2=0.45, alpha3=1-alpha1-alpha2, min.frac = 0.05, verbose = FALSE,
           lambda.accuracy=1e-4){

    if(is.null(lambdas)){
      lambdas <- betterPathCalcNew(data = data, group.index = group.index, subgroup.index = subgroup.index, alpha1 = alpha1,
                                   alpha2 = alpha2, alpha3=alpha3, min.frac = min.frac, nlam = nlam,
                                   type = "linear", lambda.accuracy = lambda.accuracy)
    }

    X <- data$x
    y <- data$y
    n <- nrow(X)
    p <- ncol(X)

    ## Setting up group lasso stuff ##

    ord.group <- order(group.index)
    group.index <- group.index[ord.group]
    groups <- unique(group.index)
    num.groups <- length(groups)

    range.group.ind <- rep(0,(num.groups+1))
    for(i in 1:num.groups){
      range.group.ind[i] <- min(which(group.index == groups[i])) - 1
    }
    range.group.ind[num.groups+1] <- ncol(X)

    group.length <- diff(range.group.ind)

    ## Setting up sub-group lasso stuff ##

    ord.subgroup <- NULL
    new.subgroup.index <- NULL
    subgroups <- NULL
    num.subgroups <- NULL
    range.subgroup.ind <- NULL
    subgroup.length <- NULL
    tmp2 <- 0

    for(k in 1:nrow(subgroup.index)){
      ## Order subgroup indices
      ord.subgroup.tmp <- order(subgroup.index[k,],na.last=NA)
      ord.subgroup <- c(ord.subgroup,ord.subgroup.tmp)

      ## Re-order subgroup index
      subgroup.index.tmp <- subgroup.index[k,ord.subgroup.tmp]
      new.subgroup.index <- c(new.subgroup.index,subgroup.index.tmp)

      ## Determine unique subgroups in each group
      subgroups.tmp <- unique(subgroup.index.tmp)
      subgroups <- c(subgroups,subgroups.tmp)

      ## Number of subgroups in each group
      num.subgroups.tmp <- length(unique(subgroup.index.tmp))
      num.subgroups <- c(num.subgroups,num.subgroups.tmp)


      ## Range.subgroup.ind
      range.subgroup.ind.tmp <- NULL

      ## To get correct range.subgroup.ind.tmp

      for(i in 1:num.subgroups.tmp){
        tmp <- min(which(subgroup.index.tmp==subgroups.tmp[i]))-1 + tmp2
        range.subgroup.ind.tmp <- c(range.subgroup.ind.tmp,tmp)
      }
      range.subgroup.ind.tmp <- c(range.subgroup.ind.tmp,length(subgroup.index.tmp)+tmp2)
      range.subgroup.ind <- c(range.subgroup.ind,range.subgroup.ind.tmp[-length(range.subgroup.ind.tmp)])

      ## length of subgroup
      subgroup.length <- c(subgroup.length, diff(range.subgroup.ind.tmp))

      ## update tmp2
      tmp2 <- length(subgroup.index.tmp)+ tmp2
    }
    range.subgroup.ind <- c(range.subgroup.ind,ncol(X))

    ## Coming up with other C++ info ##

    X <- as.matrix(X[,ord.subgroup])
    unOrd <- match(1:length(ord.subgroup),ord.subgroup)



    ## DONE SETTING UP C STUFF ##

    #alpha <- sqrt(2*log(p))/(1+sqrt(2*log(num.groups)/min(group.length)) + sqrt(2*log(p)))

    nlam = length(lambdas)
    beta.old <- rep(0,ncol(X))
    beta.is.zero.group <- rep(1,num.groups)
    beta.is.zero.subgroup <- rep(1,sum(num.subgroups))

    beta <- array(0, c(ncol(X),nlam))

    eta <- rep(0,n)

    for(k in 1:nlam){
      beta.is.zero.group <- rep(1, num.groups)
      beta.is.zero.subgroup <- rep(1,sum(num.subgroups))
      beta.old <- rep(0, ncol(X))
      eta <- rep(0,n)

      junk <- .C("CWrapper", X = as.double(as.vector(X)), y = as.double(y),
                 indexGroup = as.integer(group.index),indexSubGroup = as.integer(new.subgroup.index),
                 nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)),
                 numGroup = as.integer(num.groups),numSubGroup = as.integer(num.subgroups), totalSubGroup = as.integer(sum(num.subgroups)),
                 rangeGroupInd = as.integer(range.group.ind), rangeSubGroupInd = as.integer(range.subgroup.ind),
                 groupLen = as.integer(group.length), subGroupLen = as.integer(subgroup.length),
                 lambda1 = as.double(lambdas[k]*alpha1), lambda2 = as.double(lambdas[k]*alpha2), lambda3=as.double(lambdas[k]*alpha3),
                 beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter),
                 thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta),
                 gamma = as.double(gamma),
                 betaIsZeroGroup = as.integer(beta.is.zero.group), betaIsZeroSubGroup = as.integer(beta.is.zero.subgroup),
                 step = as.double(step), reset = as.integer(reset))

      beta.new <- junk$beta
      beta[,k] <- beta.new
      beta.is.zero.group <- junk$betaIsZeroGroup
      beta.is.zero.subgroup <- junk$betaIsZeroSubGroup
      eta <- junk$eta
      beta.old <- beta.new
      if(verbose == TRUE){
        write(paste("***Lambda", k, "***"),"")
      }
    }
    return(list(beta = beta[unOrd,], lambdas = lambdas))
  }


#######################################
# Function to find good lambda values #
#######################################


betterPathCalcNew <- function(data, group.index, subgroup.index, alpha1 = 0.45,  alpha2=0.45, alpha3=1-alpha1-alpha2,
                              min.frac = 0.05, nlam = 20, type = "linear",lambda.accuracy=1e-4){

  n <- nrow(data$x)

  if(type == "linear"){
    X <- data$x
    resp <- data$y
    n <- nrow(X)
    p <- ncol(X)
  }

  ## function which computes x_+
  threshold.function <- function(x){
    if(x<0){
      return(0)
    } else {
      return(x)
    }
  }

  ## function which computes euclidean norm
  euclidean.norm <- function(x){
    out <- sqrt(sum(x^2))
    return(out)
  }

  ## function which computes condition for which each group has it's coefficients equal to zero.
  main.condition <- function(lambda){

    group.condition.LHS <- rep(NA,nrow(subgroup.index))
    group.condition.RHS <- group.condition.LHS


    for(i in 1:nrow(subgroup.index)){
      group.condition.RHS[i] <- sum(!is.na(subgroup.index[i,]))
      subgroups    <- unique(subgroup.index[i,])
      subgroups    <- subgroups[!is.na(subgroups)]

      subgroup.condition <- rep(NA,length(subgroups))

      for(m in 1:length(subgroups)){
        ind <- subgroups[m]
        X.fit <- data.frame(X[,which(subgroup.index[i,]==ind)])
        subgroup.length <- ncol(X.fit)

        res   <- t(X.fit) %*% resp/n				# We divide by n to keep consistency with C++ program,
        # We could make CHANGES here!!! (i.e., remove n)
        soft.threshold.vector <- apply(abs(res)-alpha3*lambda,1,threshold.function)
        soft.threshold.norm   <- euclidean.norm(soft.threshold.vector)
        subgroup.condition[m] <- threshold.function(soft.threshold.norm - alpha2 * lambda * sqrt(subgroup.length))
      }

      group.condition.LHS[i] <- euclidean.norm(subgroup.condition)
    }

    out <- ( group.condition.LHS <= alpha1 * lambda * sqrt(group.condition.RHS) )
    return(out)
  }


  ## get initial range of lambda values by gradually increasing the interval

  lambda.range <- c(0,2^0)		# initial range of lambda values
  change.lambda.range <- TRUE		# should we adjust lambda range?
  j <- 0

  while(change.lambda.range){
    j <- j+1
    lambda <- lambda.range[2]
    condition.check <- main.condition(lambda)
    if(sum(condition.check)!=length(condition.check)){
      # All groups have NOT satisfied the condition, so we update lambda.range
      lambda.range[1] <- lambda
      lambda.range[2] <- 2^j
    } else {
      # All groups HAVE satisfied the condition, so we are DONE.
      change.lambda.range <- FALSE
    }
  }

  ## improve accuracy of lambda.range
  num.bisect <- 0

  while(abs(diff(lambda.range))>lambda.accuracy){
    num.bisect <- num.bisect + 1
    lambda <- sum(lambda.range)/2
    condition.check <- main.condition(lambda)
    if(sum(condition.check)!=length(condition.check)){
      # All groups have NOT satisfied the condition, so we update lambda.range
      lambda.range[1] <- lambda
    } else {
      lambda.range[2] <- lambda
    }
  }

  max.lam <- max(lambda.range)
  min.lam <- min.frac*max.lam
  lambdas <- exp(seq(log(max.lam),log(min.lam), (log(min.lam) - log(max.lam))/(nlam-1)))

  return(lambdas)
}


