#'  Distributed random projection
#' @param data is  sparse random projection matrix
#' @param K is  the number of distributed nodes.
#' @param nk is the size of subsets.
#' @param d is the dimension number.
#'        n is the sample size.
#'        p the number of variables.
#' @usage Drp(data,K, nk,d)
#' @return MSEXrp, MSEvrp, MSESrp, kopt
#' @export
#' @examples
#' K=20; nk=50; nr=10; p=8; d=5; n=K*nk;
#' data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
#' data=matrix(rpois((n-nr)*p,1),ncol=p); rexp(nr*p,1); rchisq(10000, df = 5);
#' Drp(data=data,K=K, nk=nk,d=d)
#' @importFrom stats cov
Drp=function(data,K, nk,d){
  n=nrow(data);p=ncol(data)
  th=MSEv=MSES=c(rep(1,K))
  R=matrix(rep(0,n*nk),ncol=n)
  mr=matrix(rep(0,nk*nk),ncol=nk)
  Rp=matrix(rep(0,p*d),ncol=d)
  X0=data
  for (i in 1:K){
    mr[i,]=sample(1:n,nk,replace=F);
    r=matrix(c(1:nk,sort(mr[i,])),ncol=nk,byrow=T)
    R[t(r)]=1
    X=R%*%X0
    mrp=sample(1:nk,d,replace=F);
    rp=matrix(c(mrp,1:d),ncol=d,byrow=T)
    s =max(sqrt(p),3)
    Rp= matrix(sample(x=c(sqrt(s),0,-sqrt(s)),prob=c(1/(2*s),1-(1/s),1/(2*s)),
                      replace=TRUE,size=(p*p)),nrow=p)
    #spare random projection
    Xhatr=X %*%Rp %*% t(Rp)
    th[i]=(frobenius.norm(scale(X-Xhatr)))^2/n^2
    sigmar=cov(X)
    sigmar2hat=cov(Xhatr)
    vr=eigen(sigmar)$vectors
    vrhat=eigen(sigmar2hat)$vectors
    MSEv[i]=frobenius.norm(t(vr-vrhat)%*%(vr-vrhat))/(nk);
    MSES[i]=(frobenius.norm(sigmar-sigmar2hat))^2/(n*p)
  }
  kopt=which.min(th)
  return(c(MSEXrp=min(th),MSEvrp=MSEv[kopt],MSESrp=MSES[kopt],kopt=kopt))
}

