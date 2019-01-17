#' @title Fitting a Log Cumulative Probability Model
#'
#' @description \code{lcpm} provides the maximum likelihood estimate for ordinal outcomes and a Generalized Linear Model (GLM) with the log link without the assumption of proportionality. That is, lcpm determines the MLE for log[P(y <= j)]= cut_j + X beta_j subject to [cut_{j-1} + X beta_{j-1} <= cut_j + X beta_j] and [cut_j + X beta_j <=0]. This implementation uses \code{\link{constrOptim}}  to determine the MLE and so the results account for the restricted parameter space.
#' @param formula.linear an object of class "formula": a symbolic description of the linear model to be fitted.
#' @param data dataframe containing the data in linear model.
#' @param conf.level optional confidence level (1-alpha) defaulted to 0.95.
#' @param y.order optional if y contains ordered integer categories 1:J. If y is not ordered integer 1:J then this is a vector with the ordinal values for y ranging from the lowest to largest ordinal outcome. See Examples below.
#' @param startval optional vector of the starting values.
#' @param less.than.0 optional logical for constraint cut_j <= 0 for all j=1:(J-1). Default is TRUE.
#' @param control.list optional list of controls for constrOptim
#' @param eps.outer option for constrOptim
#' @param ... Additional arguments for built in functions
#' @importFrom stats constrOptim model.frame model.matrix na.omit pchisq pnorm qnorm runif coef printCoefmat
#' @importFrom numDeriv hessian jacobian
#' @importFrom plyr count
#' @export
#' @seealso \code{\link{ppm}}
#' @return list of class "lcpm" is returned containing:
#' \item{coefficients}{vector of the estimate of cut_j and beta_j}
#' \item{se}{vector of the estimate of standard errors}
#' \item{vcov}{matrix of the inverse of the negative Hessian}
#' \item{fitted.values}{matrix of unique covariates and the corresponding estimate of the cumulative probabilities: exp(X \%*\% coefficients)}
#' \item{loglik}{numerical value of the log-likelihood at the maximum likelihood estimate}
#' \item{barrier.value}{value of mu in the log-barrier algorithm}
#' \item{outer.iterations}{value of the number of outer iterations}
#' \item{formula}{formula in the call of lcpm}
#' \item{startvalues}{vector of the starting values for constrained optimization algorithm}
#' \item{proptest}{Score test if a proportionality assumption is appropriate, includes test statistic (teststat), p-value (pval), df, and fitted proportional probability model (propmodel)}
#' @note A warning of MLE close to the boundary must be carefully considered. Data may have some structure that requires attention. Additionally, there is no imputation. Any NA results in complete row removal.
#' @references Singh, G; Fick, G.H. Ordinal Outcomes: Cumulative Probability Model with a Log Link. Manuscript in preparation.
#' @examples
#' # Example below showing the use of y.order if outcome is not integers 1:J.
#' # See examples in ppm for an additional example
#'
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
#' # LCPM estimates for non-proportional model
#' test1<-lcpm(y ~ var_a + var_b, data=testdata, y.order=c(2,3,5,10))
#' summary(test1)
#'
#' # The proportional probability model used for the score test
#' summary(test1$proptest$propmodel)

lcpm<-function(formula.linear, data,conf.level=0.95,y.order=NULL, startval=NULL, less.than.0=TRUE , control.list = NULL, eps.outer=NULL,...){

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

		if(!is.null(ncol(X))){p=ncol(X)
			} else {p=1}

		#### number of rows, columns
		Xy<-sumdat

		#nr=nrow(Xy)
		varnames<-names(Xy)[-c(1,ncol(Xy))]
		varname1<-c("start")
		for(i in 1:ncuts){
			for(j in 1:p){
			varname1<-c(varname1,paste0(varnames[j],i))
			}
		}
		varname2<-c(paste0("cut_",1:ncuts),varname1[-1],"freq")


		X1<-Xy[Xy$y==miny,]
		XJ<-Xy[Xy$y==maxy,]
		Xj<-Xy[Xy$y!=miny & Xy$y!=maxy,]

		l1<-length(X1$y)
		Xa1<-cbind(matrix(c(rep(1,l1),rep(0,l1*(ncuts-1))),ncol=ncuts,byrow=FALSE),X1[,-c(1,ncol(X1))],matrix(0,ncol=p*(ncuts-1),nrow=l1),X1[,ncol(X1)])
		Xa1<-as.matrix(Xa1)
		colnames(Xa1)<-varname2

		lJ<-length(XJ$y)
		XaJ<-cbind(matrix(c(rep(0,lJ*(ncuts-1)),rep(1,lJ)),ncol=maxy-1,byrow=FALSE),matrix(0,ncol=p*(ncuts-1),nrow=lJ),XJ[,-1])
		XaJ<-as.matrix(XaJ)
		colnames(XaJ)<-varname2

		Xj<-as.matrix(Xj)
		lj<-length(Xj[,1])
		Xaj1<-matrix(0,ncol=(p+1)*ncuts+1,nrow=lj)
		Xaj2<-matrix(0,ncol=(p+1)*ncuts+1,nrow=lj)

		for(i in 1:lj){
			cut0<-c(rep(0,ncuts))
			cut0[Xj[i,1]]<-1
			bet0<-c(rep(0,p*ncuts))
			bet0[(p*(Xj[i,1]-1)+1):(p*Xj[i,1])]<-c(Xj[i,-c(1,ncol(Xj))])
			Xaj1[i,]<-c(cut0,bet0,Xj[i,ncol(Xj)])
			cut1<-c(rep(0,ncuts))
			cut1[Xj[i,1]-1]<-1
			bet1<-c(rep(0,p*ncuts))
			bet1[(p*(Xj[i,1]-2)+1):(p*(Xj[i,1]-1))]<-c(Xj[i,-c(1,ncol(Xj))])
			Xaj2[i,]<-c(cut1,bet1,Xj[i,ncol(Xj)])
		}

			if (less.than.0==TRUE){
			extraconst=cbind(diag(rep(1,ncuts)),matrix(0,ncol=ncuts*p+1,nrow=ncuts))
			Xa<-rbind(Xa1,XaJ,Xaj1,Xaj2,-(Xaj1-Xaj2),extraconst)

			} else {
			Xa<-rbind(Xa1,XaJ,Xaj1,Xaj2,-(Xaj1-Xaj2))
			}

			Xa.f<-as.matrix(Xa[,-ncol(Xa)])

			if(is.null(startval)){
				alphastarts<-c(-1,rep(-0.1,ncuts-1))
				startval<-c(rev(cumsum(alphastarts)),rep(0,ncuts*p))
			}

			if(is.null(control.list)){control.list=list(maxit=100000,abstol=1e-12,reltol=1e-12)}
			if(is.null(eps.outer)){eps.outer=1e-10}

			if(less.than.0==TRUE){
			constr_optim1<-constrOptim(startval,lcpmMinusloglik, Xa1=Xa1, XaJ=XaJ, Xaj1=Xaj1, Xaj2=Xaj2 ,ui=-Xa.f, ci=rep(0,length(Xa.f[,1])), method="Nelder-Mead",control=control.list,outer.eps=eps.outer)
		  constr_optim<-try(constrOptim(c(constr_optim1$par[1:ncuts],rep(0,ncuts*p)),lcpmMinusloglik, Xa1=Xa1, XaJ=XaJ, Xaj1=Xaj1, Xaj2=Xaj2 ,ui=-Xa.f, ci=rep(0,length(Xa.f[,1])), method="Nelder-Mead",control=control.list,outer.eps=eps.outer))
			  if(class(constr_optim)=="try-error"){
			  constr_optim<-constr_optim1
			  }

		  } else {
				constr_optim<-constrOptim(startval,lcpmMinusloglik, Xa1=Xa1, XaJ=XaJ, Xaj1=Xaj1, Xaj2=Xaj2 ,ui=-Xa.f, ci=rep(0,length(Xa.f[,1])), method="Nelder-Mead",control=control.list,outer.eps=eps.outer)
			}

			mle<-constr_optim$par
			names(mle)<-varname2[-c(length(varname2))]
			npar<- length(mle)

			mle_alpha <- mle[1:(length(mle)-p)]

			fitted.val <- exp(rbind(Xa1,Xaj1)[,-ncol(Xa)]%*%mle)
			if(any(fitted.val>0.9999)){boundaryissue<-TRUE
			                          print("WARNING: MLE possibly on a boundary")
			} else { boundaryissue<-FALSE}
      fit.vals<- cbind(rbind(Xa1,Xaj1)[,-ncol(Xa)],fitted.val)
      colnames(fit.vals)[ncol(fit.vals)]<-"exp(X'Beta)"

			NegHess <- try(hessian(lcpmMinusloglik,constr_optim$par,Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2))
			if(class(NegHess)=="try-error"){warning("WARNING: POTENTIAL BOUNDARY ISSUE, SE ARE QUESTIONABLE")
			  testa<-matrix(NA,ncol=length(mle),nrow=length(mle))
			  SE <-c(rep(NA,length(mle)))
			} else{
              testa <- try(solve(NegHess))
              if(class(testa)!="try-error"){
                SE <-sqrt(diag(testa))
              } else {warning("WARNING: POSSIBLE SINGULAR HESSIAN MATRIX OR POTENTIAL BOUNDARY ISSUE")
                SE <-c(rep(NA,length(mle)))
              }
			}

			startval2<-as.vector(c(startval[1:ncuts],rep(0,p)))
			app0<-ppm(formula.linear,data,conf.level,y.order,startval2,less.than.0,control.list,eps.outer,...)
			app1<-ppm(formula.linear,data,conf.level,y.order,startval=c(app0$coefficients),less.than.0,control.list,eps.outer,...)
			mle_pr = as.vector(app1$coefficients)
			alpha_mle = mle[1:ncuts]


			mle_o = c(mle_pr[1:ncuts],rep(mle_pr[-c(1:ncuts)],ncuts))
			chi.stat2=t(as.vector(-jacobian(lcpmMinusloglik,mle_o,Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2)))%*%solve(hessian(lcpmMinusloglik,mle_o,Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2))%*%as.vector(-jacobian(lcpmMinusloglik,mle_o,Xa1=Xa1,XaJ=XaJ,Xaj1=Xaj1,Xaj2=Xaj2))
			p.val2 = pchisq(chi.stat2,df=((ncuts-1)*p),lower.tail=FALSE)

			out<-list(coefficients = mle, se = SE, vcov = testa, fitted.values = fit.vals,
			         loglik = -constr_optim$value,boundaryissue=boundaryissue ,convergence = constr_optim$convergence, barrier.value = constr_optim$barrier.value, outer.iterations = constr_optim$outer.iterations,
			        formula = formula.linear, startvalues = startval, proptest=list(teststat=as.numeric(chi.stat2),pval=as.numeric(p.val2),df=((ncuts-1)*p),propmodel=app1),call=match.call())
		  class(out)<-"lcpm"
		  return(out)
}
