##' Standardized Ball Imformation
##'
##' Calculate Ball Information statistic
##'
##'
##' @param x a numeric matirx.
##' @param y a numeric matirx.
##' @param R When R=0, it will return Ball information
##' @param seed Not need
##' @param weight Not need
##' @useDynLib SBISIS
##' @export
##' @examples
##'
##' n = 100
##' x = rnorm(n)
##' y = rnorm(n)
##' SBI(x, y)
##'
BI <- function(x, y, R=0, seed=2015, weight=FALSE){
	x<-as.matrix(x);y<-as.matrix(y)
	dim_x<-dim(x); dim_y<-dim(y);n<-dim_x[1]; p<-dim_x[2]; q<-dim(y)[2]
	RCT <- numeric(1)
	Dx <- numeric(n*n)
	Dy <- numeric(n*n)
	dst <- TRUE
	Dx <- .C("distance", as.double(t(x)), as.double(t(Dx)), as.integer(n),as.integer(p))
	x <- matrix(Dx[[2]],n,n)
	Dy <- .C("distance", as.double(t(y)), as.double(t(Dy)), as.integer(n),as.integer(q))
	y <- matrix(Dy[[2]],n,n)
	RCT<-.C("BI", as.double(t(x)), as.double(t(y)), as.integer(n), as.integer(p), as.integer(q), as.integer(dst), HRC=as.double(RCT), as.integer(seed), as.integer(R), as.integer(weight))
	return(RCT$HRC)
}




##' Standardized Ball Imformation
##'
##' Calculate Ball Information statistic
##'
##'
##' @param x a numeric matirx.
##' @param y a numeric matirx.
##' @export
##' @examples
##'
##' n = 100
##' x = rnorm(n)
##' y = rnorm(n)
##' SBI(x, y)
##'
SBI <- function(x,y){
  sqrt(BI(y,x,R=0)/sqrt(BI(x,x,R=0)*BI(y,y,R=0)))
}







##' Ball Information in survival
##'
##' Calculate Ball Information statistic in survival
##'
##'
##' @param x ordered covariate
##' @param t ordered survival event time
##' @param delta ordered survival event status
##' @param Sc Survfit object
##' @param n Sample size
##' @useDynLib SBISIS
##' @export
##' @examples
##'
##' library(survival)
##' data(GSE)
##' time=Y[,1]
##' delta=Y[,2]
##' Sc=survfit(Surv(time,1-delta)~1)$surv
##' xo=X[,1]
##' n=length(time)
##' #
##' ord.t = sort(time)
##' ix = order(time)
##' ord.delta = delta[ix]
##' ord.x=xo[ix]
##' SBI.surv<-SRCT(ord.x, ord.t, ord.delta, Sc, n)
##' SBI.surv
##'
SRCT <- function(x, t, delta, Sc, n){
  RCT <- numeric(1)
  RCT<-.C("SRCT", as.double(t(x)), as.double(t(t)), as.double(t(delta)),
          as.double(t(Sc)), as.integer(n), RC=as.double(RCT))
  return(RCT$RC)
}

