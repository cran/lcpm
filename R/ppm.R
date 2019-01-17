#' @title Fitting a Proportional Probability Model
#'
#'
#' @description \code{ppm} provides the maximum likelihood estimate for ordinal outcomes and a Generalized Linear Model with the log link with the assumption of proportionality. That is, ppm determines the MLE for log[P(y <= j)]= cut_j + X beta subject to [cut_{j-1} <= cut_j ] and [cut_j + X beta <=0]. This implementation uses \code{\link{constrOptim}} to determine the MLE and so the results should correctly account for the restricted parameter space. A proposed test for proportionality is included in \code{\link{lcpm}}.
#' @param formula.linear an object of class "formula": a symbolic description of the linear model to be fitted.
#' @param data dataframe containing the data in linear model.
#' @param conf.level optional confidence level (1-alpha) defaulted to 0.95.
#' @param y.order optional if y contains ordered integer categories 1:J. If y is not ordered integer 1:J then this is a vector with the ordinal values for y ranging from the lowest to largest ordinal outcome. See Examples below.
#' @param startval optional vector of the starting values.
#' @param less.than.0 optional logical for constraint cut_j <= 0 for all j=1:(J-1). Default is TRUE.
#' @param control.list optional list of controls for constrOptim.
#' @param eps.outer option for constrOptim.
#' @param ... Additional arguments for built in functions.
#' @importFrom stats constrOptim model.frame model.matrix na.omit pnorm qnorm pchisq coef printCoefmat
#' @importFrom numDeriv hessian jacobian
#' @importFrom plyr count
#' @importFrom Matrix rankMatrix
#' @export
#' @seealso \code{\link{lcpm}}
#' @return list of class "ppm" is returned containing:
#' \item{coefficients}{vector of the estimate of cut_j and beta}
#' \item{se}{vector of the estimate of standard errors}
#' \item{vcov}{matrix of the inverse of the negative Hessian}
#' \item{fitted.values}{matrix of unique covariates and the corresponding estimate of the cumulative probabilities: exp(X \%*\% coefficients)}
#' \item{loglik}{numerical value of the log-likelihood at the maximum likelihood estimate}
#' \item{barrier.value}{value of mu in the log-barrier algorithm}
#' \item{outer.iterations}{value of the number of outer iterations}
#' \item{formula}{formula in the call of ppm}
#' \item{startvalues}{vector of the starting values for constrained optimization algorithm}
#' @note A warning of MLE close to the boundary must be carefully considered. Data may have some structure that requires attention.  Additionally, there is no imputation. Any NA results in complete row removal.
#' @references Singh, G; Fick, G.H. (submitted) Ordinal Outcomes: A Cumulative Probability Model with the Log Link and the Assumption of Proportionality. Statistics in Medicine.
#' @examples
#' # 2 examples below showing the use of y.order if outcome are not integers 1:J.
#'
#' # Example 1:
#' 
#' var_a <- c(rep(0,60),rep(1,60))
#' var_b <- c(rep(0,90),rep(1,30))
#' y1<-c(rep(2,5),rep(3,10),rep(5,5),rep(10,10),
#' rep(2,5),rep(3,10),rep(5,10),rep(10,5),
#' rep(2,10),rep(3,5),rep(5,5),rep(10,10),
#' rep(2,10),rep(3,5),rep(5,10),rep(10,5))
#'
#' testdata<-data.frame(y=y1,var_a=var_a,var_b=var_b)
#' 
#' # PPM estimates for proportional model
#' test1<-ppm( y ~ var_a + var_b, data=testdata, y.order=c(2,3,5,10))
#' summary(test1)
#'
#' # Example 2:
#'
#' y2<-c(rep("a",5),rep("b",10),rep("c",5),rep("d",10),
#' rep("a",5),rep("b",10),rep("c",10),rep("d",5),
#' rep("a",10),rep("b",5),rep("c",5),rep("d",10),
#' rep("a",10),rep("b",5),rep("c",10),rep("d",5))
#' testdata2<-data.frame(y=y2,var_a=var_a,var_b=var_b)
#' test2<-ppm(y~var_a + var_b , data=testdata2, y.order=c("a","b","c","d"))
#' summary(test2)

ppm<-function (formula.linear, data,conf.level=0.95, y.order = NULL, startval = NULL, less.than.0=TRUE, control.list = NULL ,eps.outer=NULL,...){

	if(!is.data.frame(data)){
		print("Make data into a dataframe")
		stop()
		}

	if(conf.level >= 1 | conf.level<=0){
		print("conf.level must be greater than 0 but less than 1")
		stop()
		}

	Y<-model.frame(formula.linear, data)[, 1]
	X<- model.matrix(formula.linear, data)

	if(!is.null(ncol(X))){p=ncol(X)-1
		} else {p=1}

	if(!is.null(y.order)) {
		newY <- match(Y,y.order)
		mydata<-na.omit(as.data.frame(cbind(newY,X[,-1])))
		colnames(mydata)<-c("y",colnames(X)[-1])
	}
	else{mydata<-na.omit(as.data.frame(cbind(Y,X[,-1])))
			colnames(mydata)<-c("y",colnames(X)[-1])}

  sumdat<-count(mydata)


						miny=1
						maxy=max(sumdat$y)
						ncuts=length(unique(mydata$y))-1
						y<-sumdat[,1]
						X<-sumdat[,-c(1,ncol(sumdat))]

						Xy<-sumdat

						X1<-Xy[Xy$y==miny,]
						XJ<-Xy[Xy$y==maxy,]
						Xj<-Xy[Xy$y!=miny & Xy$y!=maxy,]

						concavity.XJ<-cbind(rep(1,nrow(XJ)),XJ[,-c(1,ncol(XJ))])
						if(p==1){strictconcavity=TRUE
						} else if(rankMatrix(concavity.XJ)==ncol(concavity.XJ)){strictconcavity=TRUE
						} else{strictconcavity=FALSE
						}

						l1<-length(X1$y)
						Xa1<-cbind(matrix(c(rep(1,l1),rep(0,l1*(maxy-2))),ncol=maxy-1,byrow=FALSE),X1[,-1])
						Xa1<-as.matrix(Xa1)

						lJ<-length(XJ$y)
						XaJ<-cbind(matrix(c(rep(0,lJ*(maxy-2)),rep(1,lJ)),ncol=maxy-1,byrow=FALSE),XJ[,-1])
						XaJ<-as.matrix(XaJ)

						lj<-length(Xj$y)
						Xaj1<-cbind(matrix(0,ncol=maxy-1,nrow=lj),Xj)
						Xaj1[cbind(seq(1,lj,1),Xaj1$y)]<-1
						Xaj1<-as.matrix(Xaj1[,-maxy])

						Xaj2<-cbind(matrix(0,ncol=maxy-1,nrow=lj),Xj)
						Xaj2[cbind(seq(1,lj,1),Xaj2$y-1)]<-1
						Xaj2<-as.matrix(Xaj2[,-maxy])

						if(less.than.0==TRUE){
							extraconst=cbind(diag(rep(1,ncuts)),matrix(0,ncol=p+1,nrow=ncuts))
							Xa<-rbind(Xa1,XaJ,Xaj1,Xaj2,-(Xaj1-Xaj2),extraconst)
							} else {
									Xa<-rbind(Xa1,XaJ,Xaj1,Xaj2,-(Xaj1-Xaj2))
							}



						Xa.f<-as.matrix(Xa[,-ncol(Xa)])

						if(is.null(startval)){
						alphastarts<-c(-1,rep(-0.1,ncuts-1))
						startval<-c(rev(cumsum(alphastarts)),rep(0,p))
						}

						if(is.null(control.list)){control.list=list(maxit=100000,abstol=1e-12,reltol=1e-12)}
						if(is.null(eps.outer)){eps.outer=1e-10}

						if(less.than.0==TRUE){
						constr_optim1<-constrOptim(startval,lcpmMinusloglik, Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2,ui=-Xa.f, ci=rep(0,length(Xa.f[,1])), method="Nelder-Mead",control=control.list,outer.eps=eps.outer)
						constr_optim<-constrOptim(c(constr_optim1$par[1:ncuts],rep(0,p)),lcpmMinusloglik, Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2,ui=-Xa.f, ci=rep(0,length(Xa.f[,1])), method="Nelder-Mead",control=control.list,outer.eps=eps.outer)
						} else {
							constr_optim<-constrOptim(startval,lcpmMinusloglik, Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2,ui=-Xa.f, ci=rep(0,length(Xa.f[,1])), method="Nelder-Mead",control=control.list,outer.eps=eps.outer)
						}



			mle<-constr_optim$par
			mle_alpha <- mle[1:(length(mle)-p)]
			names(mle)<-c(paste("cut_",1:length(mle_alpha),sep=""),names(XJ)[-c(1,ncol(XJ))])
			npar<- length(mle)

			fitted.val <- exp(rbind(Xa1,Xaj1)[,-ncol(Xa)]%*%mle)
			if(any(fitted.val>0.9999)){boundaryissue<-TRUE
			print("WARNING: MLE possibly on a boundary")
			} else { boundaryissue<-FALSE}
			fit.vals<- cbind(rbind(Xa1,Xaj1)[,-ncol(Xa)],fitted.val)
			colnames(fit.vals)[ncol(fit.vals)]<-"exp(X'Beta)"


			NegHess<-try(hessian(lcpmMinusloglik,constr_optim$par,Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2))
			if(class(NegHess)=="try-error"){warning("WARNING: POTENTIAL BOUNDARY ISSUE, SE ARE QUESTIONABLE")
		    testa<-matrix(NA,ncol=length(mle),nrow=length(mle))
		    SE <-c(rep(NA,length(mle)))
			} else{

      			testa <- try(solve(NegHess))
      			if(class(testa)!="try-error"){
      			  SE <-sqrt(diag(testa))
      			}
      			else {warning("WARNING: POSSIBLE SINGULAR HESSIAN MATRIX OR POTENTIAL BOUNDARY ISSUE")
      			  SE <-c(rep(NA,length(mle)))}
			}

			out<-list(coefficients = mle, se = SE, vcov = testa, fitted.values = fit.vals,Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2,
			         loglik = -constr_optim$value,boundaryissue=boundaryissue,strictconcavity=strictconcavity , convergence = constr_optim$convergence, barrier.value = constr_optim$barrier.value, outer.iterations = constr_optim$outer.iterations,
			         formula = formula.linear, startvalues = startval,call=match.call())

      class(out)<-"ppm"
      return(out)
}
