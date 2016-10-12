\name{SBI.sis}
\alias{SBI.sis}
\title{Standardized Ball Information Sure Independence Screening}
\description{Standardized Ball Information Sure Independence Screening}
\usage{
SBI.sis(X, Y, candidate = c("large"), method = c("SBI-SIS", "SBI-IISIS-lm", "SBI-IISIS-gam"), parms = list(d1 = 5, d2 = 5, df = 3))
}
\arguments{
  \item{X}{a numeric matirx.}
  \item{Y}{a numeric matirx.}
  \item{candidate}{size of candidate set, large mean sample size, small mean sample size divide logarithm sample size}
  \item{method}{Method for screening procedure, include "SBI-SIS", "SBI-IISIS-lm", "SBI-IISIS-gam"."SBI-SIS" means Sure Independence Screening Procedure, "SBI-IISIS-lm" means Iterative "SBI-SIS" using linear regression, "SBI-IISIS-gam" means Iterative "SBI-SIS" using generalized additive models}
  \item{parms}{Parameters used to control the Iterative "SBI-SIS". d1 is the size of initial set, d2 is the variable set size added in each iteration}
}


\examples{
n<-30
p<-100
error<-rnorm(n,0,1)
creat.sigma1 <- function(rho,p) {
  Sigma<-matrix(0,p,p)
  for(i in 1:p)
    for(j in 1:p)
      Sigma[i,j]=rho^abs(i-j)
  return(Sigma)
}
sigma = creat.sigma1(0.5,p)
ev <- eigen(sigma, symmetric = TRUE)
if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
  warning("sigma is numerically not positive definite")
}
Sigma <- ev$vectors \%*\% diag(sqrt(ev$values), length(ev$values)) \%*\% t(ev$vectors)
vsubset=c(1,2,3,4)
X<-matrix(rnorm(n*p),n,p)\%*\%Sigma
beta<-matrix(c(5,5,5,-15),4,1)
Y<-as.vector(cbind(X[,1],X[,2],X[,3],X[,4])\%*\%beta)+error
re=SBI.sis(X = X,Y = Y,candidate = "large",method = "SBI-IISIS-gam",parms = list(d1=5,d2=5,df=4))
re$candidate.set
}

