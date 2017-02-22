\name{SBISIS-package}
\alias{SBISIS-package}
\alias{SBISIS}
\docType{package}
\title{Standardized Ball Information Sure Independence Screening}
\description{This R package is about standardized ball information sure independence screening and its extensions}
\details{
\tabular{ll}{
Package: \tab SBISIS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-10-13\cr
License: \tab GPL-2\cr
}
The package contains two main function:

\code{SBI.sis} - Standardized ball information sure independence screening(SBI-SIS) and interative SBI-SIS

\code{SBI.sis.surv} - Extension of SBI-SIS in survival data

}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
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
