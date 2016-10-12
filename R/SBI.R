BI.test <- function(x, y, R=0, seed=2015, weight=FALSE){
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

SBI.sis<-function(X,Y,candidate=c("large"),method=c("SBI-SIS","SBI-IISIS-lm","SBI-IISIS-gam"),parms=list(d1=5,d2=5,df=3))
{
  n=dim(X)[1];p=dim(X)[2]
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
  d1=parms$d1
  d2=parms$d2
  df=parms$df
  if(candidate=="small"){
    final_d=d_logn
  } else {
    final_d=d
  }

  rcory_result=c()

  #first round
  rcory_result=lapply(as.list(1:length(ids)),function(id,x=X,y=Y){
    xo=x[,id]
    sqrt(BI.test(y,xo,R=0)/sqrt(BI.test(xo,xo,R=0)*BI.test(y,y,R=0)))
  })
  # rnorm()
  rcory_result=unlist(rcory_result)
  max_ids=order(rcory_result,decreasing=T)

  method_sub=paste0(head(unlist(strsplit(method,"-")),2),collapse='-')
  method_sub1=tail(unlist(strsplit(method,"-")),1)
  if(method_sub=="SBI-IISIS"){
    chooseids=max_ids[1:d1]
    Xhavepickout=ids[chooseids]
    Xlastpickout=ids[chooseids]  # flag the pick out variables, used to remove the effect of selected variables

    ids=ids[-chooseids]
    #iteration round
    if(method_sub1=='lm'){
      while(length(Xhavepickout)<final_d)
      {
        # lm fit for X
        Xnew=lm(X[,ids]~X[,Xhavepickout])$resid
        # lm fit for Y
        Y=lm(Y~X[,Xlastpickout])$resid

        # SBI-screening
        rcory_result=lapply(as.list(1:length(ids)),function(id,x=Xnew,y=Y){
          xo=x[,id]
          sqrt(BI.test(y,xo,R=0)/sqrt(BI.test(xo,xo,R=0)*BI.test(y,y,R=0)))
        })
        rcory_result=unlist(rcory_result)
        max_ids=order(rcory_result,decreasing=T)

        # prepare for next iteration
        chooseids=max_ids[1:d2]
        Xhavepickout=c(Xhavepickout,ids[chooseids])
        Xlastpickout=ids[chooseids]
        ids=ids[-chooseids]
      }
    }
    if(method_sub1=='gam'){
      while(length(Xhavepickout)<final_d)
      {
        # gam fit for X
        lastpickout_formula=paste0(' + s(',colnames(X)[Xlastpickout],collapse = paste0(",df = ",df,")"))
        lastpickout_formula=paste0(lastpickout_formula,paste0(",df = ",df,")"),collapse = "")
        lastpickout_dat=X[,Xlastpickout]
        Xnew=sapply(ids,function(x){
          formula_one=paste0(colnames(X)[x],"~",lastpickout_formula)
          formula_one=as.formula(formula_one)
          dat=as.data.frame(cbind(X[,x],lastpickout_dat))
          colnames(dat)[1]=colnames(X)[x]
          # dat=as.data.frame(dat)
          # colnames(dat)=paste0("X",c(x,Xhavepickout))
          gam(formula_one,data = dat)$residuals
        })

        # gam fit for Y
        dat=data.frame(Y,lastpickout_dat)
        names(dat)[1]=c("Y")
        formula_Y=as.formula(paste("Y~",lastpickout_formula))
        Y=gam(formula = formula_Y,data = dat)$residuals

        # SBI-screening
        rcory_result=lapply(as.list(1:length(ids)),function(id,x=Xnew,y=Y){
          xo=x[,id]
          sqrt(BI.test(y,xo,R=0)/sqrt(BI.test(xo,xo,R=0)*BI.test(y,y,R=0)))
        })
        rcory_result=unlist(rcory_result)
        max_ids=order(rcory_result,decreasing=T)

        # prepare for next iteration
        chooseids=max_ids[1:d2] #
        Xhavepickout=c(Xhavepickout,ids[chooseids])
        Xlastpickout=ids[chooseids]
        ids=ids[-chooseids]
      }
    }
  } else {
    chooseids=max_ids[1:final_d]
    Xhavepickout=ids[chooseids]
  }
  return(list(candidate.set=Xhavepickout,
              candidate.data=list(X=X[,Xhavepickout],Y=Y),
              candidate.size=length(Xhavepickout)))
}

SRCT <- function(x, t, delta, Sc, n){
  RCT <- numeric(1)
  RCT<-.C("SRCT", as.double(t(x)), as.double(t(t)), as.double(t(delta)),
          as.double(t(Sc)), as.integer(n), RC=as.double(RCT))
  return(RCT$RC)
}

SBI.sis.surv<-function(X,Y,candidate=c("large")){

  n=dim(X)[1];p=dim(X)[2]
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

  # SBI Screening(survival)
  fitc = survfit(Surv(time,1-delta)~1)
  Sc = fitc$surv
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

