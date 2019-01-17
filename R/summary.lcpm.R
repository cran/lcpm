#' @export
#'
summary.lcpm<-function(object, ...){

  se <- object$se
  zval <- coef(object)/se
  TAB <- cbind(Estimate = coef(object), StdErr = se, z.value = zval,
               p.value = 2 * pnorm(-abs(zval), mean = 0, sd = 1))


  res <- list(call = object$call, coefficients = TAB, loglik = object$loglik,boundaryissue=object$boundaryissue,
              proptest = object$proptest$teststat, propval = object$proptest$pval, propdf = object$proptest$df)
  class(res) <- "summary.lcpm"
  res
}
