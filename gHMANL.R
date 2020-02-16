#############################################################################################################################
################     Gene-based high-dimensional nonlinear mediation analysis (gHMA-NL)
#############################################################################################################################
#' Gene-based high-dimensional nonlinear mediation analysis (gHMA-NL)
#'
#' @param X The exposure vector. Each component corresponds to a sample.
#' @param Y The outcome vector. Each component corresponds to a sample.
#' @param M The observed methylation matrix. Rows represent samples and columns are CpG sites.
#' @param COV.XM The matrix of covariates dataset for testing the association X ~ M. Default = NULL.
#' @param COV.MY The matrix of covariates dataset for testing the association Y ~ M. Using covariates should be careful. If the cavariables are not specified, the covariates for Y ~ M are the same with that of M ~ X.
#' @param family Either 'gaussian' or 'binomial', depending on the data type of outcome (Y). 
#' @param thr the proportion of variance explained in the KPCA to determine the number of KPCs.
#'
#' @return A data.frame containing the following items. 
#' \item{pcom_a}{a gene-level p-value of the α arm effect (from exposure to mediators).}
#' \item{pcom_b}{a gene-level p-value of the β arm effect (from mediators to outcome adjusting for exposure).}
#' \item{pcom}{a gene-level p-value of the total mediation effect.}
#' 
library(KRLS)
library(Matrix)
library(rms)
source("function_FisherCombPval.R")

gHMANL <- function(X, Y, M, COV.XM = NULL, COV.MY = COV.XM,
                   family = c("gaussian", "binomial"), thr = 0.80)
{
  #############################################################################################################################
  ################################   Stage 1  KPCA -- gaussian kernel          #################################
  #############################################################################################################################
  ker  <- gausskernel(X=M,sigma=ncol(M)) # Gaussian kernel
  ker  <- as.matrix(forceSymmetric(ker))
  ker0 <- ker
  nn   <- nrow(M)
  diag(ker0) <- rep(0,nn) 
  J <- matrix(1,nn,nn)
  ker.cen <- ker-J%*%ker0/(nn-1)-ker0%*%J/(nn-1)+J%*%ker0%*%J/nn/(nn-1)
  kc <- ker.cen
  
  res      <- eigen(kc/nn, symmetric=TRUE)
  features <- which(res$values > 0)
  values   <- res$values[features]
  pcv      <- t(t(res$vectors[, features])/sqrt(res$values[features]))
  cumVar   <- cumsum(values)/sum(values)
  n_kpc    <- min(which(cumVar >= thr)[1], ncol(M))
  kpcs     <- kc %*% pcv[, 1:n_kpc]
  
  kpcs   <- as.matrix(kpcs)
  rownames(kpcs) <- rownames(M)
  colnames(kpcs) <- paste0("pc_", 1:ncol(kpcs))
  
  #############################################################################################################################
  ################################   Stage 2  mediation analysis          #################################
  #############################################################################################################################
  M <- kpcs
  n <- nrow(M)
  p <- ncol(M)
  if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2 * n/log(n)) 
  d <- min(p, d) # if d > p select all mediators
  #############################################################################################################################
  ################################           Step-1 Assess the α arm effect          #################################
  #############################################################################################################################
  ps=pl=pc=rep(0, ncol(M))
  for (k in 1 : ncol(M))
  {
    if(is.null(COV.XM)){
      pl[k] <- summary(lm(M[,k]~X))$coefficients[2,4] 
      pc[k] <- summary(lm(M[,k]~X))$coefficients[2,1]
    }else{
      pl[k] <- summary(lm(M[,k]~X + COV.XM))$coefficients[2,4] 
      pc[k] <- summary(lm(M[,k]~X + COV.XM))$coefficients[2,1]
    }
    if(is.null(COV.MY)){
      ps[k] <- summary(glm(Y~M[,k] + X, family = family))$coefficients[2,4] 
    }else{
      ps[k] <- summary(glm(Y~M[,k] + X + COV.MY, family = family))$coefficients[2,4] 
    }
  }
  
  ID      <- which(ps<=sort(ps)[d])
  M_SIS   <- M[, ID]
  XM      <- cbind(M_SIS, X)
  ID_test <- ID
  
  alpha_est <- t(cbind(pc, pl))
  colnames(alpha_est) <- colnames(M)
  alpha_est <- alpha_est[,ID_test, drop = FALSE]
  alpha_hat <- as.numeric(alpha_est[1, ])         #  the estimator for alpha
  pval_a    <- alpha_est[2, ]
  
  ####### combine the p-values of alpha by conducting Fisher combination test proposed by Yang et al.
  if(length(pval_a)==1){
    p.com_a  <- pval_a
  } else{
    pheno     <- M[, pmatch(names(pval_a), colnames(M))]
    pval      <- t(as.matrix(pval_a))
    p.com_a   <- FisherCombPval(pval, pheno)
  }
  #############################################################################################################################
  ################################           Step-2 Assess the β arm effect          #################################
  #############################################################################################################################
  if (is.null(COV.MY)) {
    YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
  }
  res <- summary(glm(Y ~ .,family = family, data = YMX))$coefficients
  
  beta_est <- t(cbind(res[2:(length(ID_test) + 1), 1], res[2:(length(ID_test) + 1), 4]))
  beta_hat <- beta_est[1,]
  P_beta   <- beta_est[2,]
  pval_b   <- P_beta
  
  ####### conduct likelihood ratio test (LRT) 
  if (is.null(COV.MY)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV.MY)
  }
  
  if(length(pval_b)==1){
    p.com_b  <- pval_b
  } else{
    library(rms)
    fit.0 <- glm(Y ~ ., family = family, data = YX)
    fit.1 <- glm(Y ~ ., family = family, data = YMX)
    lrt <- lrtest(fit.0, fit.1)
    p.com_b <- as.numeric(lrt[[1]][3])
  }
  
  #############################################################################################################################
  ################################           Step-3 Assess the total mediation effect          #################################
  #############################################################################################################################
  ####### using the joint intersection-union test
  p.com <- max(p.com_a, p.com_b)
  
  result <- data.frame(pcom_a=p.com_a, pcom_b=p.com_b, pcom=p.com)
  
  return(result)
}
