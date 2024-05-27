\name{FPC}
\alias{FPC}
\title{FPC}
\usage{
FPC(data,K,nk)}
\description{
FPC is a technology that applies FPC A to distributed computing environments.
}
\arguments{
\item{data}{is  a data set matrix.}
\item{K}{is  the desired target rank.}
\item{nk}{is the size of subsets.}

}
\value{
MSEv1,MSEv2,MSEvopt,MSESopt1,MSESopt2,MSESopt,MSEShat,MSESba,MSESw}
\examples{
K=20; nk=500; p=8; n=10000;m=50
data=matrix(c(rnorm((n-m)*p,0,1),rpois(m*p,100)),ncol=p)
FPC(data=data,K=K,nk=nk)}


