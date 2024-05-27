#'  Distributed random svd
#' @param data is  sparse random projection matrix
#' @param K  is the number of distributed nodes.
#' @param nk is the size of subsets.
#' @param m  the dimension of variables
#' @param q number of additional power iterations.
#' @param k the desired target rank.
#' @return MSE of Xrsvd,vrsvd,rsvd and kopt.
#' @examples
#'install.packages("matrixcalc")#install.packages("rsvd")
#'library(matrixcalc)
#'library(rsvd)
#'K=20; nk=50; nr=10; p=8; m=5; q=5;k=4;n=K*nk;
#' data=X=matrix(rexp(n*p,0.8),ncol=p)
#' Drsvd(data=data,K=K,nk=nk,m=m,q=q,k=k)
Drsvd=function(data,K, nk,m,q,k){
  n=nrow(data);p=ncol(data); q=q;
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
    Vk=rsvd(X,m,p,q)$v; Zk=X%*%Vk
    #Vk=rpca(X,m,p,q)$rotation[,1:k]; Zk=X%*%Vk;#rpca
    # Vk=rpca(rrpca(X)$L,m,p,q)$rotation[,1:k]; Zk=X%*%Vk; #rrpca
    #Vk=rsvd(rqb(X,m,p,q)$B,m,p,q)$v[,1:k]; Zk=X%*%Vk;#rqb
    Xhatr=Zk%*%t(Vk)#rsvd
    th[i]=(frobenius.norm(scale(X-Xhatr)))^2/n^2
    sigmar=cov(X)
    sigmar2hat=cov(Xhatr)
    vr=eigen(sigmar)$vectors
    vrhat=eigen(sigmar2hat)$vectors
    MSEv[i]=frobenius.norm(t(vr-vrhat)%*%(vr-vrhat))/(nk);
    MSES[i]=(frobenius.norm(sigmar-sigmar2hat))^2/(n*p)
  }
  kopt=which.min(th)
  return(c(MSEXrsvd=min(th),MSEvrsvd=MSEv[kopt],MSESrsvd=MSES[kopt],kopt=kopt))
}

