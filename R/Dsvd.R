#' @param data a real input matrix (or data frame) to be decomposed.
#' @param K the number of blocks into which variable X is divided.
#' @param nk The number of each blocks.
#' @param k  the desired target rank.
#' @return MSE of Xs,vsvd,Ssvd and kopt.
#' @export
#' @examples
#'install.packages("matrixcalc")
#'library(matrixcalc)
#' K=20; nk=50; nr=10; p=8; k=4; n=K*nk;
#' data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
#' Dsvd(data=data,K=K, nk=nk,k=k)
Dsvd=function(data,K, nk,k){
  n=nrow(data);p=ncol(data)
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
    v=svd(X)$v
    vk=v[,1:k]
    Xhats=X%*%vk%*%t(vk)
    S0=cov(scale(X0))
    Msvd=cov(scale(Xhats))
    th[i]=(frobenius.norm(scale(X-Xhats)))^2/n^2
    sigmas=cov(X)
    sigmas2hat=cov(Xhats)
    vs=eigen(sigmas)$vectors
    vshat=eigen(sigmas2hat)$vectors
    MSEv[i]=frobenius.norm(t(vs-vshat)%*%(vs-vshat))/(nk);
    MSES[i]=(frobenius.norm(sigmas-sigmas2hat))^2/(n*p)
  }
  kopt=which.min(th)
  return(c(MSEXs=min(th),MSEvsvd=MSEv[kopt],MSESsvd=MSES[kopt],kopt=kopt))
}

