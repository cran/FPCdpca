\name{Dpca}
\alias{Dpca}
\title{Distributed PCA}
\usage{
Dpca(data,K, nk)
}
\description{
Distributed PCA is a technology that applies PCA to distributed computing environments.
}
\arguments{
\item{data}{is the n random vectors constitute the data matrix.}
\item{K}{is an index subset/sub-vector specifying.}
\item{nk}{is the size of subsets.}
}
\value{
MSEXp, MSEvp, MSESp, kopt
}
\examples{
K=20; nk=50; nr=10; p=8;n=K*nk;d=6
data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
Dpca(data,K,nk)
}

