\name{SRCT}
\alias{SRCT}
\title{Ball Information of survival time response}
\description{Computes Ball Information statistic of survival time response}
\usage{SRCT(x, t, delta, Sc, n)}

\arguments{
  \item{x}{ordered covariate}
  \item{t}{ordered survival event time}
  \item{delta}{ordered survival event status}
  \item{Sc}{Survfit object}
  \item{n}{Sample size}
}

\examples{
library(survival)
data(GSE)
time=Y[,1]
delta=Y[,2]
Sc=survfit(Surv(time,1-delta)~1)$surv
xo=X[,1]
n=length(time)
#
ord.t = sort(time)
ix = order(time)
ord.delta = delta[ix]
ord.x=xo[ix]
SBI.surv<-SRCT(ord.x, ord.t, ord.delta, Sc, n)
SBI.surv
}

