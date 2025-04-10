#' Distributed Random SVD
#'
#' @param data A numeric matrix or data frame.
#' @param K Number of distributed nodes.
#' @param nk Size of each subset.
#' @param m Target dimension for random projection.
#' @param q Number of power iterations.
#' @param k Desired rank.
#' @return A vector containing MSE values and optimal subset index.
#' @export
#' @examples
#' library(rsvd)
#' library(matrixcalc)
#' K <- 20
#' nk <- 50
#' p <- 8
#' m <- 5
#' q <- 5
#' k <- 4
#' n <- K * nk
#' data <- matrix(rexp(n * p, 0.8), ncol = p)
#' Drsvd(data = data, K = K, nk = nk, m = m, q = q, k = k)
#' @importFrom stats cov
Drsvd <- function(data, K, nk, m, q, k) {
  n <- nrow(data)
  p <- ncol(data)
  X0 <- data
  th <- MSEv <- MSES <- rep(1, K)
  R <- matrix(0, nrow = n, ncol = n)
  mr <- matrix(0, nrow = nk, ncol = nk)
  Rm <- matrix(0, nrow = nk, ncol = K)
  
  for (i in 1:K) {
    mr[i, ] <- sample(1:n, nk, replace = FALSE)
    r <- matrix(c(1:nk, sort(mr[i, ])), ncol = nk, byrow = TRUE)
    Rm[, i] <- r[2, ]
    R[t(r)] <- 1
    X <- R %*% X0
    Vk <- rsvd::rsvd(X, m, p, q)$v  # 显式调用 rsvd
    Zk <- X %*% Vk
    Xhatr <- Zk %*% t(Vk)
    th[i] <- (matrixcalc::frobenius.norm(scale(X - Xhatr)))^2 / n^2
    sigmar <- stats::cov(X)
    sigmar2hat <- stats::cov(Xhatr)
    vr <- eigen(sigmar)$vectors
    vrhat <- eigen(sigmar2hat)$vectors
    MSEv[i] <- matrixcalc::frobenius.norm(t(vr - vrhat) %*% (vr - vrhat)) / nk
    MSES[i] <- (matrixcalc::frobenius.norm(sigmar - sigmar2hat))^2 / (n * p)
  }
  kopt <- which.min(th)
  return(c(
    MSEXrsvd = min(th),
    MSEvrsvd = MSEv[kopt],
    MSESrsvd = MSES[kopt],
    kopt = kopt
  ))
}