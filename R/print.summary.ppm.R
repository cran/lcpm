#' @export

print.summary.ppm<-function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n")
  if(x$strictconcavity){cat("loglikelihood: strictly concave in the interior")
    }else{cat("WARNING: loglikelihood concave in the interior")}
  cat("\n")
  cat("loglikelihood: ",x$loglik)
  cat("\n")
  if(x$boundaryissue){cat("WARNING: MLE possibly on a boundary")
  cat("\n")}

}
