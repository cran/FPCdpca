#'  Decentralized PCA
#'
#' @param data is  a sparse random projection matrix
#' @param K is  the desired target rank.
#'        nk is the size of subsets.
#'        eps
#'        nit.max
#'        n is the sample size.
#'        p the number of variables.
#'
#' @return MSEXrp, MSEvrp, MSESrp, kopt
#' @export
#'
#' @examples
#' K=20; nk=50; nr=10; p=8;  n=K*nk;d=5
#' data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
#' set.seed(1234)
#' eps=10^(-1);nit.max=1000
#' Depca(data=data,K=K, nk=nk,eps=eps,nit.max=nit.max)
#' TXde=TSde=c(rep(0,5))
#' for (j in 1:5){
#'  depca=Depca(data,K, nk,d,k,eps,nit.max)
#'  TXde[j]=as.numeric(depca)[1]
#'  TSde[j]=as.numeric(depca)[2]}
#' mean(TXde)
#' mean(TSde)

Depca=function(data,K,nk, eps,nit.max){
  d=6
  n=nrow(data);p=ncol(data)
  niter=0;
  X0=data
  th=MSEv=MSES=c(rep(1,K))
  R=matrix(rep(0,n*nk),ncol=n)
  mr=matrix(rep(0,nk*nk),ncol=nk)
  Rm=matrix(rep(0,nk*K),ncol=K)
  for (i in 1:K){
    mr[i,]=sample(1:n,nk,replace=F);
    r=matrix(c(1:nk,sort(mr[i,])),ncol=nk,byrow=T)
    Rm[,i]=r[2,]
    R[t(r)]=1
    X=R%*%X0
    A=t(X)%*%X
    E=cov(X)
    eV=eigen(E)$vectors
    S=W=eV[,1:d]
    diff=10
    while(diff > eps & niter < nit.max){ Wold=W
    Sold=S
    W=qr.Q(qr(S))
    S=S+A%*%(W-Wold)
    diff=min(norm((S-Sold),type="2"))
    niter=niter+1
    }
    XdePCA=X%*%W%*%t(W)
    th[i]=(frobenius.norm(scale(X-XdePCA)))^2/n^2
    sigmad=cov(X)
    sigmad2hat=cov(XdePCA)
    vd=eigen(sigmad)$vectors
    vdhat=eigen(sigmad2hat)$vectors
    MSEv[i]=frobenius.norm(t(vd-vdhat)%*%(vd-vdhat))/(nk);
    MSES[i]=(frobenius.norm(sigmad-sigmad2hat))^2/(n*p)
  }
  kopt=which.min(th)
  return(c(MSEXrp=min(th),MSEvrp=MSEv[kopt],MSESrp=MSES[kopt],kopt=kopt))
}
