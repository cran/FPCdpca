\name{Depca}
\alias{Depca}
\title{Decentralized PCA}
\usage{
Depca(data,K,nk, eps,nit.max)
}
\description{
Decentralized PCA
is a technology that applies decentralized PCA  to distributed computing environments.
}
\arguments{
\item{data}{is  sparse random projection matrix.}
\item{K}{is  the desired target rank.}
\item{nk}{is the size of subsets.}
\item{eps}{is the noise.}
\item{nit.max}{is the repeat times.}
}
\value{
MSEXrp,MSEvrp, MSESrp, kopt
}
\examples{
K=20; nk=50; nr=10; p=8; k=4; n=K*nk;d=6
data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
set.seed(1234)
eps=10^(-1);nit.max=1000
TXde=TSde=c(rep(0,5))
for (j in 1:5){
  depca=Depca(data=data,K=K, nk=nk,eps=eps,nit.max=nit.max)
  TXde[j]=as.numeric(depca)[1]
  TSde[j]=as.numeric(depca)[2]
}
mean(TXde)
mean(TSde)
}


