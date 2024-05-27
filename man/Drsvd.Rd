\name{Drsvd}
\alias{Drsvd}
\title{ Distributed random svd}
\usage{
Drsvd(data,K, nk,m,q,k)
}
\description{
Distributed random svd is a technology that applies random SVD to distributed computing environments.
}
\arguments{
\item{data}{ sparse random projection matrix.}

\item{K}{the number of distributed nodes.}

\item{nk}{the size of subsets.}

\item{m}{the dimension of variables.}

\item{q}{number of additional power iterations.}

\item{k}{the desired target rank.}
}

\value{
 \item{MSEXrsvd}{The MSE value of Xrsvd}
 \item{MSEvrsvd}{The MSE value of vrsvd}
 \item{MSESrsvd}{The MSE value of Srsvd}
 \item{kopt}{The size of optimal subset}
}
\examples{
K=20; nk=50; nr=10; p=8; m=5; q=5;k=4;n=K*nk;
data=X=matrix(rexp(n*p,0.8),ncol=p)
#data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
#data=X=matrix(rpois((n-nr)*p,1),ncol=p); rexp(nr*p,1); rchisq(10000, df = 5);
#data=X=matrix(rexp(n*p,0.8),ncol=p)
Drsvd(data=data,K=K,nk=nk,m=m,q=q,k=k)
}

