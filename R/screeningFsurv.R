##' Standardized Ball Information Sure Independence Screening
##'
##' Standardized Ball Information Sure Independence Screening
##'
##'
##' @param X a numeric matirx.
##' @param Y a numeric matirx or Surv object
##' @param candidate size of candidate set, large mean sample size, small mean
##' sample size divide logarithm sample size
##' @param standize Whether standize the explantory matrix?(Default TRUE)
##' @import survival
##' @export
##' @examples
##'
##' data(GSE)
##' re=SBI.sis.surv(X = X,Y = Y,candidate = "large")
##' re$candidate.set
##'
SBI.sis.surv<-function(X, Y, candidate=c("large"), standize = TRUE){

  n=dim(X)[1]; p=dim(X)[2]
  n = as.numeric(n); p = as.numeric(p)
  ids=1:p
  Y=as.matrix(Y)
  X=as.matrix(X)
  colnames(X)=paste0("X",1:p)
  colnames(Y)=paste0("Y",1:ncol(Y))
  if(any(apply(Y,2,anyNA))) {stop("NA appear in matrix Y")}
  if(any(apply(X,2,anyNA))) {stop("NA appear in matrix X")}

  # decide candicate size
  d_logn=round(n/log(n))
  d=n
  if(candidate=="small"){
    final_d=d_logn
  } else {
    final_d=d
  }

  # prepare for screening
  time=Y[,1]
  delta=Y[,2]
  ord.t = sort(time)
  ix = order(time)
  ord.delta = delta[ix]
  xo=X[ix,]
  if(standize) {
    xo = apply(xo, 2, scale)
  }

  # SBI Screening(survival)
  fitc = survfit(Surv(time,1-delta)~1)
  Sc = fitc$surv
  if(length(unique(ord.t)) != n) {
    rep_num = as.data.frame(table(ord.t))[, "Freq"]
    Sc = mapply(function(x, y) {
      rep(x, y)
    }, Sc, rep_num, SIMPLIFY = FALSE)
    Sc= unlist(Sc)
  }
  # num_non_zero = n - sum()
  rcory_result=apply(xo,2,function(x){
    SRCT(x = x,t = ord.t,delta = ord.delta,Sc = Sc,n = n)
  })
  rcory_result=unlist(rcory_result)
  max_ids=order(rcory_result,decreasing=T)
  chooseids=max_ids[1:final_d]
  Xhavepickout=ids[chooseids]

  return(list(candidate.set=Xhavepickout,
              candidate.data=list(X=X[,Xhavepickout],Y=Y),
              candidate.size=length(Xhavepickout)))
}
