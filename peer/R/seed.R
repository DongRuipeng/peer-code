seed.thresholding <- function(S, u0, theta, niter = 400) {
  u.new <- u0
  for (i in c(1:niter)) {
    u.old <- u.new
    u <- S %*% u.old
    ind <- abs(u) < theta
    u[ind] <- 0
    if (norm(u, "2") == 0) break
    u.new <- u / norm(u, "2")
    if (norm(u.old - u.new, "2") < 1e-6) break
  }
  u <- u.new
  return(u)
}

seed.opt <- function(X, Y, theta, nrank, rho = 1e-2) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  P <- t(X) %*% X / n + rho * diag(p)
  P.inv <- solve(P)
  U <- matrix(0, p, nrank);
  V <- matrix(0, q, nrank);
  d <- rep(0, nrank)
  Y.temp <- Y
  for (k in c(1:nrank)) {
    R <- t(Y.temp) %*% X / n
    Q <- t(R) %*% R / q
    S <- P.inv %*% Q
    eig <- eigen(S, symmetric = TRUE)
    u0 <- eig$vectors[, 1]
    u <- seed.thresholding(S, u0, theta)
    u <- u / norm(X %*% u, "2")
    v <- t(Y.temp) %*% X %*% u / (norm(X %*% u, "2") ^ 2)
    Y.temp <- Y.temp - X %*% u %*% t(v)
    d[k] <- norm(v, "2");
    U[, k] <- u;
    V[, k] <- v / d[k]
  }
  C = U %*% diag(d) %*% t(V)
  mse <- norm(Y - X %*% C, "F") ^ 2 / (n*q)
  GIC <- sqrt(n) * log(mse) + sqrt(log(p * q)) * log(log(n)) * sum(U != 0)
  temp <- list(C = C, U = U, V = V, d = d, GIC = GIC)
  return(temp)
}


#' Sequential estimation with eigen-decompsition (SEED) 
#' @param X predictors matrix
#' @param Y responses matrix 
#' @param nrank rank of estimator 
#' @param ntheta number of thresholding parameters for tuning parameter 
#' @return 
#'  \item{C}{estimator for coefficient}
#'  \item{U}{estimator for left singular vectors}
#'  \item{V}{estimator for right singular vectors}
#'  \item{d}{estimator for singular values}
#' @export
seed <- function(X, Y, nrank, ntheta = 100) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  theta.range <- c(log(0.5), log(norm(t(Y) %*% X / n, "M")))
  theta <- c(exp(seq(theta.range[1], theta.range[2], length = ntheta)))
  U.stor <- array(0, c(p, nrank, ntheta)) 
  V.stor <- array(0, c(q, nrank, ntheta))
  d.stor <- matrix(0, nrank, ntheta)
  gic.stor <- rep(0, ntheta)
  C.stor <- array(0, c(p, q, ntheta))
  for(i in c(1:ntheta)){
    fit <- seed.opt(X, Y, theta[i], nrank)
    C.stor[,,i] <- fit$C 
    U.stor[,,i] <- fit$U 
    V.stor[,,i] <- fit$V 
    d.stor[,i] <- fit$d 
    gic.stor[i] <- fit$GIC 
  }
  index <- sort(gic.stor, index.return = TRUE)$ix[1]
  C <- C.stor[,,index]; U <- U.stor[,,index]; V <- V.stor[,,index]
  d <- d.stor[,index] 
  return(list(C = C, U = U, V = V, d = d)) 
}