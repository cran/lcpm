#' @export
#'
summary.ppm<-function(object, ...){

  se <- object$se
  zval <- coef(object)/se
  TAB <- cbind(Estimate = coef(object), StdErr = se, z.value = zval,
               p.value = 2 * pnorm(-abs(zval), mean = 0, sd = 1))


  res <- list(call = object$call, coefficients = TAB, loglik = object$loglik,boundaryissue=object$boundaryissue,strictconcavity=object$strictconcavity)
  class(res) <- "summary.ppm"
  res
}
