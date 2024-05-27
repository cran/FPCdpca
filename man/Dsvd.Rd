\name{Dsvd}
\alias{Dsvd}
\title{Distributed svd}
\usage{
Dsvd(data,K, nk,k)
}
\description{
Distributed svd is a technology that applies SVD to distributed computing environments.
}
\arguments{
\item{data}{A independent variable.}
\item{K}{the number of distributed nodes.}
\item{nk}{the number of each blocks.}
\item{k}{the desired target rank.}
}
\value{
 \item{MSEXs}{the MSE of Xs}
 \item{MSEvsvd}{the MSE of vsvd}
 \item{MSESsvd}{the MSE of Ssvd}
 \item{kopt}{the size of optimal subset}
}
\examples{
#install.packages("matrixcalc")
library(matrixcalc)
K=20; nk=50; nr=10; p=8; k=4; n=K*nk;
data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
Dsvd(data=data,K=K, nk=nk,k=k)
}


