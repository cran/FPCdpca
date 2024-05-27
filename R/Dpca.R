#'  The distributed PCA
#'
#' data is the n random vectors constitute the data matrix.
#'  K is an index subset/sub-vector specifying.
#'        nk is the size of subsets.
#'        n is the sample size.
#'        p the number of variables.
#'
#' @return MSEXp, MSEvp, MSESp, kopt
#' @export
#'
#' @examples
#' K=20; nk=50; nr=10; p=8;  n=K*nk;d=6
#' data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
#' Dpca(data=data,K=K, nk=nk)

Dpca=function(data,K, nk){
  n=nrow(data);p=ncol(data)
  X0=data
  d=6
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
    E=cov(X)
    eV=eigen(E)$vectors
    evd=eV[1:d,]
    Xhatp=X%*%t(evd) %*%evd
    th[i]=(frobenius.norm(scale(X-Xhatp)))^2/n^2
    sigmap=cov(X)
    sigmap2hat=cov(Xhatp)
    vp=eigen(sigmap)$vectors
    vphat=eigen(sigmap2hat)$vectors
    MSEv[i]=frobenius.norm(t(vp-vphat)%*%(vp-vphat))/(nk);
    MSES[i]=(frobenius.norm(sigmap-sigmap2hat))^2/(n*p)
  }
  kopt=which.min(th)
  return(c(MSEXp=min(th),MSEvp=MSEv[kopt],MSESp=MSES[kopt],kopt=kopt))
}
