

#####################################
## Data generation for Group lasso ##
#####################################
## N 	    : sample size
## p.group : number of covariates in each group
## p.subgroup : matrix of number of covariates in each subgroup (each row corresponds to a group)
## beta.coef : matrix of coefficient vectors (rows: groups, columns: coefficient values)
## sigma   : model error covariance
## index   : vector indicating which covariates are in which groups

#' Data generation for group lasso
#'
#' @param N sample size
#' @param p.group number of covariates in each subgroup
#' @param beta.coef matrix of coefficient vectors (rows: groups, columns: coefficient values)
#' @param sigma model error variance
#'
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix bdiag
#' @import stats
#' @export
data.group <- function(N,p.group,beta.coef,sigma){

  ########################################
  ## Form beta vector ##
  ########################################

  beta <- as.vector(t(beta.coef))	## Make into a vector
  beta <- beta[!is.na(beta)]		## Remove NA values

  ##################################################################
  ## Make p.subgroup into a vector ##
  ##################################################################

  p.group.vector <- p.group

  ######################################
  ## Form covariates ##
  ######################################
  ## X \sim Normal(mu,Sigma)

  no.subgroup <- length(p.group.vector)
  Sigma.list <- vector("list",no.subgroup)
  for(j in 1:no.subgroup){
    Sigma.list[[j]] <- 0.7 * matrix(1,nrow=p.group.vector[j],
                                    ncol=p.group.vector[j]) +
      0.3 * diag(1,nrow=p.group.vector[j])
  }

  ## Form Sigma covariance matrix
  X.Sigma <- as.matrix(Matrix::bdiag(Sigma.list))

  ## Form mean of normal distribution
  X.mu <- rep(0,length(beta))

  ## Form covariates
  X <- MASS::mvrnorm(N,X.mu,X.Sigma)
  colnames(X) <- paste("X_",seq(1,sum(p.group)),sep="")


  ################################################
  ## Form response vector ##
  ################################################

  y <- X %*% beta + stats::rnorm(N,mean=0,sd=sigma)


  list(y=y,X=X)
}
