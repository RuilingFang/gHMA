#############################################################################################################################
################     Gene-based high-dimensional linear mediation analysis (gHMA-L)
#############################################################################################################################
#' Gene-based high-dimensional linear mediation analysis (gHMA-L)
#'
#' @param X The exposure vector. Each component corresponds to a sample.
#' @param Y The outcome vector. Each component corresponds to a sample.
#' @param M The observed methylation matrix. Rows represent samples and columns are CpG sites.
#' @param COV.XM The matrix of covariates dataset for testing the association X ~ M. Default = NULL.
#' @param COV.MY The matrix of covariates dataset for testing the association Y ~ M. Using covariates should be careful. If the cavariables are not specified, the covariates for Y ~ M are the same with that of M ~ X.
#' @param penalty The penalty to be applied to the model. Either 'MCP', 'SCAD', or 'lasso'.
#' @param family Either 'gaussian' or 'binomial', depending on the data type of outcome (Y). 
#' @param topN An integer can be used to set the number of top markers by the method of sure independent screening. Default = NULL. 
#'             If topN is NULL, it will be either ceiling(2*n/log(n)) if family = 'gaussian', or ceiling(n/(2*log(n))) if family = 'binomial', where n is the sample size.
#'             If the sample size is greater than topN (pre-specified or calculated), all markers will be harbored in the test. 
#'
#' @return A data.frame containing the following items. 
#' \item{pcom_a}{a gene-level p-value of the α arm effect (from exposure to mediators).}
#' \item{pcom_b}{a gene-level p-value of the β arm effect (from mediators to outcome adjusting for exposure).}
#' \item{pcom}{a gene-level p-value of the total mediation effect.}
#' 

library(ncvreg)
library(rms)
source("function_FisherCombPval.R")
gHMAL <- function(X, Y, M, COV.XM = NULL, COV.MY = COV.XM, 
                  penalty = c("MCP", "SCAD", "lasso"), 
                  family = c("gaussian", "binomial"), topN = NULL)
{
  n <- nrow(M)
  p <- ncol(M)
  if (is.null(topN)) {
    if (family == "binomial") d <- ceiling(n/(2*log(n))) else d <- ceiling(2 * n/log(n)) 
  } else {
    d <- topN  # the number of top mediators that associated with exposure (X)
  }
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
  
  ID <- which(ps<=sort(ps)[d])
  M_SIS <- M[, ID]
  M_SIS <- as.matrix(M_SIS)
  XM    <- cbind(M_SIS, X)
  
  if (is.null(COV.MY)) {
    fit <- ncvreg(XM, Y, family = family, penalty = penalty,
                  penalty.factor = c(rep(1, ncol(M_SIS)), 0))
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    XM_COV <- cbind(XM, COV.MY)
    fit <- ncvreg(XM_COV, Y, family = family, penalty = penalty, 
                  penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV.MY))))
  }
  
  lam <- fit$lambda[which.min(BIC(fit))]
  
  est <- coef(fit, lambda = lam)[2:(d+1)]
  ID_test <- ID[which(est!= 0)]
  if(length(ID_test) == 0) { ID_test <- ID }
  
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
  ####### fitting the 3rd equation in model (2) #########################
  if (is.null(COV.MY)) {
    YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
  }
  res <- summary(glm(Y ~ .,family = family, data = YMX))$coefficients
  
  beta_est <- t(cbind(res[2:(length(ID_test) + 1), 1], res[2:(length(ID_test) + 1), 4]))
  beta_hat <- beta_est[1, ]
  P_beta   <- beta_est[2, ]
  pval_b   <- P_beta
  
  ####### conduct likelihood ratio test (LRT) ######################
  if (is.null(COV.MY)) {
    YX <- data.frame(Y = Y, X = X)
  } else {
    YX <- data.frame(Y = Y, X = X, COV.MY)
  }
  if(length(pval_b)==1){
    p.com_b  <- pval_b
  } else{
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