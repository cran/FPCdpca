#' Distributed Principal Component Analysis (DPCA)
#'
#' Performs distributed PCA on a data matrix partitioned into subsets.
#'
#' @param data A numeric matrix or data frame containing the data, where rows are observations and columns are variables.
#' @param K Integer, the number of subsets to partition the data into.
#' @param nk Integer, the size of each subset (number of rows per subset).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{MSEXp}: Minimum squared reconstruction error.
#'   \item \code{MSEvp}: MSE of eigenvectors.
#'   \item \code{MSESp}: MSE of covariance matrix.
#'   \item \code{kopt}: Optimal subset index.
#' }
#'
#' @details
#' The function splits the input data matrix into \code{K} subsets of size \code{nk} each. 
#' The parameters \code{n} (number of rows) and \code{p} (number of columns) are automatically 
#' derived from the input data matrix as \code{n = nrow(data)} and \code{p = ncol(data)}.
#'
#' @examples
#' K <- 20
#' nk <- 50
#' nr <- 10
#' p <- 8
#' n <- K * nk
#' d <- 6
#' data <- matrix(c(rnorm((n - nr) * p, 0, 1), rpois(nr * p, 100)), ncol = p)
#' Dpca(data = data, K = K, nk = nk)
#' @importFrom stats cov
#'
#' @export

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
