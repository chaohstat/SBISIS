\name{SBI.sis.surv}
\alias{SBI.sis.surv}
\title{Ball Information Sure Independence Screening for survival time response}
\description{Computes Ball Information statistics of survival time response for sure independence screening }
\usage{
SBI.sis.surv(X, Y, candidate = c("large"))
}
\arguments{
  \item{X}{a numeric matirx.}
  \item{Y}{a numeric matirx or Surv object}
  \item{candidate}{size of candidate set, large mean sample size, small mean sample size divide logarithm sample size}
}

\examples{
data(GSE)
re=SBI.sis.surv(X = X,Y = Y,candidate = "large")
re$candidate.set
}

