#' FPC
#' @param data is a data set
#' @param K is an index subset/sub-vectoris
#' @param nk is an index subset/sub-vectoris for each block
#' @usage FPC(data,K,nk)
#' @return MSEv1,MSEv2,MSEvopt,MSESopt1,MSESopt2,MSESopt,MSEShat,MSESba,MSESw
#' @export
#' @examples
#' library(matrixcalc)
#' K=20; nk=500; p=8; n=10000;m=50
#' data=matrix(c(rnorm((n-m)*p,0,1),rpois(m*p,100)),ncol=p)
#' FPC(data=data,K=K,nk=nk)
#' @importFrom stats var
FPC=function(data,K,nk){
  n=nrow(data);p=ncol(data)
  X0=data
  Tr=tv=th=fv=hv=c(rep(1,K))
  mu=w=muni=c(rep(0,K))
  R=matrix(rep(0,n*nk),ncol=n)
  Io=matrix(rep(0,nk*K),ncol=nk);
  mr=matrix(rep(0,nk*nk),ncol=nk)
  Rm=matrix(rep(0,nk*K),ncol=K)
  S1=matrix(rep(0,p*p),ncol=p)
  M4=matrix(rep(0,p*p),ncol=p)
  S0=cov(scale(X0))
  S0hat=(t(scale(X0))%*%scale(X0))/n
  for (i in 1:K){
    oldS1=S1
    oldSw1=M4
    mr[i,]=sample(1:n,nk,replace=F);
    r=matrix(c(1:nk,sort(mr[i,])),ncol=nk,byrow=T)
    Rm[,i]=r[2,]
    R[t(r)]=1
    X=R%*%X0
    M1=cov(X)
    eigenM1=eigen(M1)
    eigenValue1=eigenM1$values
    eigenVector1=eigenM1$vectors
    orderValue1=order(eigenValue1,decreasing=T)
    v0=eigenVector1[,orderValue1]
    M2=(t(scale(X))%*%scale(X))/nk
    eigenM2=eigen(M2)
    eigenValue2=eigenM2$values
    eigenVector2=eigenM2$vectors
    orderValue2=order(eigenValue2, decreasing=T)
    v0hat=eigenVector2[,orderValue2]
    for (k in 1:K){
      mu[k]=sum(diag(var(M2)))
      muni[k]=solve(mu[k])
    }
    for (k in 1:K){
      w[k]=muni[k]/sum(muni)
    }
    th[i]=frobenius.norm(t(v0)%*%v0hat)
    fv[i]=frobenius.norm(v0hat%*%t(v0hat)-v0%*%t(v0))
    sigmaw=w[k]*M2
    M4=sigmaw+oldSw1
    S1=M2+oldS1
  }
  M3=S1/K
  opt1=Rm[,which.min(fv)]
  opt2=Rm[,which.max(th)]
  opt3=intersect(opt1,opt2) ; opt3
  X1=X0[opt1,]
  sigma1=cov(X1)
  v1=eigen(sigma1)$vectors
  sigma1hat=(t(X1-mean(X1))%*%(X1-(mean(X1))))/(nrow(X1)-1)
  v1hat=eigen(sigma1hat)$vectors
  MSEv1=frobenius.norm(t(v1-v1hat)%*%(v1-v1hat))/(nk)
  X2=X0[opt2,]
  sigma2=cov(X2)
  v2=eigen(sigma2)$vectors
  sigma2hat=(t(X2-mean(X2))%*%(X2-(mean(X2))))/(nrow(X2)-1)
  v2hat=eigen(sigma2hat)$vectors
  MSEv2=frobenius.norm(t(v2-v2hat)%*%(v2-v2hat))/(nk)
  X3=X0[opt3,]
  sigma3=cov(X3)
  v3=eigen(sigma3)$vectors
  sigma3hat=(t(X3-mean(X3))%*%(X3-(mean(X3))))/(nrow(X3)-1)
  v3hat=eigen(sigma3hat)$vectors
  MSEvopt=frobenius.norm(t(v3-v3hat)%*%(v3-v3hat))/(nk);
  ku1=sum (scale(X1)^4)/(n*p)
  ku2=sum (scale(X2)^4)/(n*p)
  ku3=sum (scale(X3)^4)/(n*p)
  MSEShat=(frobenius.norm(S0-S0hat))^2/(n*p)
  MSESba=(frobenius.norm(S0-M3))^2/(n*p)
  MSESw=(frobenius.norm(S0-M4))^2/(n*p)
  MSESopt=(frobenius.norm(sigma3-sigma3hat))^2/(n*p)
  MSESopt2=(frobenius.norm(sigma2-sigma2hat))^2/(n*p)
  MSESopt1=(frobenius.norm(sigma1-sigma1hat))^2/(n*p)
  return(c(MSEv1=MSEv1, MSEv2=MSEv2, MSEvopt=MSEvopt,       MSESopt1=MSESopt1,MSESopt2=MSESopt2,MSESopt=MSESopt,MSEShat=MSEShat,MSESba=MSESba,MSESw=MSESw))
}

