#' @export

print.summary.lcpm<-function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  cat("loglikelihood: ",x$loglik)
  cat("\n")
  if(x$boundaryissue){cat("WARNING: MLE possibly on a boundary")
    cat("\n")}
  cat("Score Test of Proportionality: ")
  cat("Chisq test statistic ", x$proptest, " with ", x$propdf, " degrees of freedom and p-value ", x$propval)
  cat("\n")

}
