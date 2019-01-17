#' @title The Negative Log-Likelihood Function for a Log Cumulative Probability Model
#'
#' @description \code{lcpmMinusloglik} provides the negative of the log-likelihood function for a Generalized Linear Model with a log link and ordinal outcomes to be minimized in functions \code{\link{lcpm}} and \code{\link{ppm}}.
#' @param betapar a vector of values.
#' @param Xa1 matrix of covariates for all subjects with the lowest ordinal outcome value 1.
#' @param XaJ matrix of covariates for all subjects with the largest ordinal outcome value J.
#' @param Xaj1 matrix of covariates for all subjects with the ordinal outcomes with value 1 < j < J.
#' @param Xaj2 matrix of covariates for all subjects with the ordinal outcome with value 1 < j < J but lagged by 1.
#' 
#' @return value of the negative log-likelihood evaluated at betapar 
#' @export

lcpmMinusloglik<-function(betapar,Xa1,XaJ,Xaj1,Xaj2){

  loglik.1<-t(Xa1[,ncol(Xa1)])%*%(Xa1[,-ncol(Xa1)]%*%betapar)
  loglik.J<-t(XaJ[,ncol(XaJ)])%*%log(1-exp(XaJ[,-ncol(XaJ)]%*%betapar))
  loglik.j<-t(Xaj1[,ncol(Xaj1)])%*%log(exp(Xaj1[,-ncol(Xaj1)]%*%betapar)-exp(Xaj2[,-ncol(Xaj2)]%*%betapar))

  loglik.i<- as.numeric(loglik.1)+as.numeric(loglik.J)+as.numeric(loglik.j)
  minusloglik<--loglik.i
  minusloglik
}
