\name{Drpca}
\alias{Drpca}
\title{Distributed random PCA}
\usage{
Drpca(data,K, nk,d)
}
\description{
Distributed random PCA is a technology that applies random PCA to distributed computing environments.
}
\arguments{
\item{data}{is sparse random projection matrix.}
\item{K}{is  the number of distributed nodes.}
\item{nk}{is the size of subsets.}
\item{d}{is the dimension number.}
}
\value{
MSEXrp, MSEvrp, kSopt, kxopt
}
\examples{
K=20; nk=50; nr=50; p=8;d=5; n=K*nk;
data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
Drpca(data,K, nk,d)
}

