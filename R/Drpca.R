#'  Distributed random PCA
#'
#' @param data is  sparse random projection matrix
#' @param K is  the number of distributed nodes.
#' @param nk is the size of subsets.
#' @param d is the dimension number.
#'        n is the sample size.
#'        p the number of variables.
#' @usage Drpca(data,K, nk,d)
#' @return MSEXrp, MSEvrp, kSopt, kxopt
#' @export
#' @examples
#' K=20; nk=50; nr=50; p=8;d=5; n=K*nk;
#' data=matrix(c(rnorm((n-nr)*p,0,1),rpois(nr*p,100)),ncol=p)
#' @importFrom stats cov
Drpca=function(data,K, nk,d){
  n=nrow(data);p=ncol(data)
  th=fv=c(rep(1,K))
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
    #Rp[rp]=1
    s =max(sqrt(p),3)
    Rp= matrix(sample(x=c(sqrt(s),0,-sqrt(s)),prob=c(1/(2*s),1-(1/s),1/(2*s)),
                      replace=TRUE,size=(p*p)),nrow=p) #spare random projection
    Xhatr=X %*%Rp %*% t(Rp)
    th[i]=(frobenius.norm(scale(X-Xhatr)))^2/n^2
    S0=cov(scale(X0))
    Mrpca=cov(scale(Xhatr))
    #fv[i]=(frobenius.norm(scale(S0)-scale(Mrpca)))^2/p^2
    fv[i]=(frobenius.norm(S0-Mrpca))^2/(n*p)
  }
  return(c(MSEXrp=min(th),MSESrp=min(fv,na.rm = TRUE),kxopt=which.min(th),
           kSopt=which.min(fv)))
}
