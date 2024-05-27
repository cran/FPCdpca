\name{Drp}
\alias{Drp}
\title{Distributed random projection}
\usage{
Drp(data,K, nk,d)
}
\description{
Distributed random projection is a technology that applies random projection to distributed computing environments.
}
\arguments{
\item{data}{is sparse random projection matrix.}
\item{K}{is the number of distributed nodes.}
\item{nk}{is the size of subsets.}
\item{d}{is the dimension number.}
}
\value{
MSEXrp,MSEvrp, MSESrp, kopt
}
\examples{
K=20; nk=50; nr=10; p=8; d=5; n=K*nk;
data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
data=matrix(rpois((n-nr)*p,1),ncol=p); rexp(nr*p,1); rchisq(10000, df = 5);
Drp(data=data,K=K, nk=nk,d=d)
}
