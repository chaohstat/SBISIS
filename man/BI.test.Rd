\name{SBISIS}
\alias{BI.test}

\title{Ball Imformation}

\description{Calculate Ball Information statistic}

\usage{
BI.test(x, y, R = 0, seed = 2015, weight=FALSE)
}

\arguments{
  \item{x}{a numeric matirx.}
  \item{y}{a numeric matirx.}
  \item{R}{the number of replication, when R equals to 0, the function return the statistic of Ball Imformation.}
  \item{seed}{the random seed.}
}

\examples{
n = 100
x = rnorm(n)
y = rnorm(n)
BI.test(x, y, R=0)
}
