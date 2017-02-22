##' Standardized Ball Information Sure Independence Screening
##'
##' Standardized Ball Information Sure Independence Screening
##'
##'
##' @param X a numeric matirx.
##' @param Y a numeric matirx.
##' @param candidate size of candidate set, large mean sample size, small mean
##' sample size divide logarithm sample size
##' @param method Method for screening procedure, include "SBI-SIS",
##' "SBI-IISIS-lm", "SBI-IISIS-gam"."SBI-SIS" means Sure Independence Screening
##' Procedure, "SBI-IISIS-lm" means Iterative "SBI-SIS" using linear
##' regression, "SBI-IISIS-gam" means Iterative "SBI-SIS" using generalized
##' additive models
##' @param parms Parameters used to control the Iterative "SBI-SIS". d1 is the
##' size of initial set, d2 is the variable set size added in each iteration
##' @import stats
##' @import gam
##' @import splines
##' @export
##' @examples
##'
##' n<-30
##' p<-100
##' error<-rnorm(n,0,1)
##' creat.sigma1 <- function(rho,p) {
##'   Sigma<-matrix(0,p,p)
##'   for(i in 1:p)
##'     for(j in 1:p)
##'       Sigma[i,j]=rho^abs(i-j)
##'   return(Sigma)
##' }
##' sigma = creat.sigma1(0.5,p)
##' ev <- eigen(sigma, symmetric = TRUE)
##' if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
##'   warning("sigma is numerically not positive definite")
##' }
##' Sigma <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
##' vsubset=c(1,2,3,4)
##' X<-matrix(rnorm(n*p),n,p)%*%Sigma
##' beta<-matrix(c(5,5,5,-15),4,1)
##' Y<-as.vector(cbind(X[,1],X[,2],X[,3],X[,4])%*%beta)+error
##' re=SBI.sis(X = X,Y = Y,candidate = "large",method = "SBI-IISIS-gam",parms = list(d1=5,d2=5,df=4))
##' re$candidate.set
##'
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
    SBI(y,xo)
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
          SBI(y,xo)
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
          SBI(y,xo)
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
