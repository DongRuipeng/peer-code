threshold.rank <- function(A, nrank) {
  A.svd <- svd(A, nu = nrank, nv = nrank)
  if(nrank == 1) {
    temp <- A.svd$d * A.svd$u %*% t(A.svd$v)
  } else {
    temp <- A.svd$u %*% diag(A.svd$d[1:nrank]) %*% t(A.svd$v)
  }
  return(temp)
}

completion <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

get.gic <- function(X, y, Beta) {
  n <- nrow(X)
  p <- ncol(X)
  nbeta <- ncol(Beta)
  gic <- rep(0, nbeta)
  for (i in c(1:nbeta)) {
    mse <- norm(y - X %*% Beta[, i], "2") ^ 2 / n
    gic[i] <- log(mse) + log(p) * log(log(n)) * sum(Beta[, i] != 0) / n
  }
  return(gic)
}

#' Estimate the true rank of C* 
#' @param Y responses matrix 
#' @param erank initial estimated rank  
#' @return 
#'  \item{\code{r}}{refined rank estimation}
#' @export 
est.rank <- function(Y, erank) {
  n <- nrow(Y)
  q <- ncol(Y)
  miss.set <- is.na(Y)
  num_miss <- sum(miss.set)

  if (num_miss == 0) {
    tilde_d <- svd(Y, nu = erank, nv = erank)$d[1:erank]
    tilde_d <- tilde_d / sqrt(n * q)
    delta <- tilde_d[1:(erank - 1)] - tilde_d[2:erank]
    tau_n <- log(log(n)) / log(n) 
    index <- c(1:(erank - 1))
    r <- max(index[delta > tau_n])
  } else {
    Y <- apply(Y, 2, completion);
    tY.new <- Y
    for (iteration in c(1:100)) {
      tY.old <- tY.new
      tY.new <- threshold.rank(Y, erank)
      if (norm(tY.old - tY.new, "F") / norm(tY.old, "F") < 1e-6) break
      Y[miss.set] <- tY.new[miss.set]
    }
    tilde_d <- svd(tY.new, nu = erank, nv = erank)$d[1:erank]
    tilde_d <- tilde_d / sqrt(n * q)
    delta <- tilde_d[1:(erank - 1)] - tilde_d[2:erank]
    m <- n * q - num_miss
    tau_n <- log(log(n)) / log(n) 
    index <- c(1:(erank - 1))
    if (length(index[delta > tau_n]) == 0) {
      r <- 1
    } else {
      r <- max(index[delta > tau_n])
    }
  }
  return(r)
}

#' Parallel estimation with extracted layers and regularization
#' @param X predictors matrix
#' @param Y responses matrix 
#' @param nrank rank of estimator 
#' @param penalty penalty type (l1, l0, default l1)
#' @param eps tolerance parameter 
#' @return 
#'  \item{C}{estimator for coefficient}
#'  \item{U}{estimator for left singular vectors}
#'  \item{V}{estimator for right singular vectors}
#'  \item{d}{estimator for singular values}
#'  \item{time}{elapsedtime in parallel version}
#' @import glmnet BeSS stats
#' @export 
peer <- function(X, Y, nrank = 3, penalty = "l1", eps = 1e-6) {
  n <- nrow(Y);
  p <- ncol(X);
  q <- ncol(Y)
  t_start <- proc.time()
  # completion step 
  miss.set <- is.na(Y)
  if (sum(miss.set) != 0) {
    Y <- apply(Y, 2, completion);
    tY.new <- Y
    for (iteration in c(1:100)) {
      tY.old <- tY.new
      tY.new <- threshold.rank(Y, nrank = nrank)
      if (norm(tY.old - tY.new, "F") / norm(tY.old, "F") < eps) break
      Y[miss.set] <- tY.new[miss.set]
    }
    Y <- tY.new
  }
  # divide step 
  layers <- svd(Y / sqrt(n), nu = nrank, nv = nrank)
  Z <- sqrt(n) * layers$u
  V <- layers$v
  d <- layers$d
  t_common <- proc.time() - t_start
  U <- matrix(0, p, nrank)
  # do in parallel 
  t = rep(0, nrank)
  for (k in c(1:nrank)) {
    if (penalty == "l1") {
      start_time <- proc.time()
      # fit <- cv.glmnet(, intercept = FALSE, standardize = FALSE) 
      fit <- glmnet(X, Z[, k], standardize = FALSE, intercept = FALSE)
      tt <- proc.time() - start_time
      Beta <- as.matrix(coef(fit))[-1,]
      gic <- get.gic(X, Z[, k], Beta)
      ind <- sort(gic, index.return = TRUE)$ix[1]
      U[, k] <- as.vector(Beta[, ind])
    } else if (penalty == "l0") {
      start_time <- proc.time()
      fit <- bess(X, Z[, k], method = "sequential")
      tt <- proc.time() - start_time
      ind <- sort(fit$GIC, index.return = TRUE)$ix[1]
      U[, k] <- as.vector(coef(fit)[, ind])[-1]
    } else {
      stop("Not support this penalty!")
    }
    t[k] <- tt[1]
  }
  if (nrank == 1) {
    C = d * U %*% t(V)
  }else{
    C = U %*% diag(d[1:nrank]) %*% t(V)
  }
  time <- max(t) + t_common[1]
  return(list(C = C, U = U, V = V, d = d[1:nrank], time = time))
}

#' generate data and coefficient 
#' @param n number of observations
#' @param p dimension of predictor vector
#' @param q dimension of reponse vector 
#' @param nrank rank of true coefficient 
#' @param s sparsity of left singular vectors
#' @param snr singnal to noise rate 
#' @param miss missing rate of response matrix
#' @return 
#'  \item{X}{predictor matrix}
#'  \item{Y}{response matrix}
#'  \item{d}{singular values}
#'  \item{U}{left singular vectors}
#'  \item{V}{right singular vectors}
#'  \item{C}{coefficient matrix}
#' @export 
#' @import stats
sim.setup <- function(n = 100, p = 100, q = 100, nrank = 3, s = 4, snr = 0.5, miss = 0) {
  U <- matrix(0, p, nrank)
  V <- matrix(0, q, nrank)
  for (k in c(1:nrank)) {
    U[, k] <- c(rep(0, (k - 1) * s), sample(c(1, -1), s, replace = TRUE), rep(0, p - k * s))
    V[, k] <- sample(c(1, -1), q, replace = TRUE) * runif(q, 0.3, 1)
    U[, k] <- U[, k] / norm(U[, k], "2")
    V[, k] <- V[, k] / norm(V[, k], "2")
  }
  d <- sort(5 + 5 * c(1:nrank), decreasing = TRUE)
  Xsigma <- 0.5 ^ abs(outer(1:p, 1:p, FUN = "-"))
  data <- gen.data(d, U, V, n, snr, Xsigma)
  X <- data$X
  Y <- data$Y
  if (miss != 0) {
    t.ind <- sample.int(n * q, size = miss * n * q)
    y <- as.vector(Y)
    y[t.ind] <- NA
    Y <- matrix(y, n, q)
  }

  temp <- list(X = X, Y = Y, d = d, U = U, V = V, C = U %*% diag(d) %*% t(V))
  return(temp)
}

#' generate data based on the given coefficients
#' @param d singular values
#' @param U left singular vectors
#' @param V right singular vectors
#' @param n sample size
#' @param snr singnal to noise rate 
#' @param Xsigma covariance matrix of predictors 
#' @return 
#'  \item{X}{predictor matrix}
#'  \item{Y}{response matrix}
#' @export
#' @import MASS stats
gen.data <- function(d, U, V, n, snr, Xsigma) {
  D <- diag(d)
  # finding basis along more number of columns of data vector 
  basis.vec <- function(x) {
    if (diff(dim(x)) < 0) x <- t(x)
    qd <- qr(x)
    k <- qr.Q(qd) %*% qr.R(qd)[, 1:qd$rank]
    k[abs(k) < 1e-6] <- 0
    b.ind <- vector()
    for (i in 1:qd$rank)
      b.ind[i] <- which(apply(x, 2, function(x, y) sum(abs(x - y)), k[, i]) < 1e-6)[1]
    return(list(ind = b.ind, vec = x[, b.ind]))
  }

  p <- nrow(U);
  q <- nrow(V);
  nrank <- ncol(U)

  U.t <- diag(max(dim(U)))
  U.t <- U.t[, - basis.vec(U)$ind]
  P <- cbind(U, U.t)
  UtXsUt <- t(U.t) %*% Xsigma %*% U.t
  UtXsU <- t(U.t) %*% Xsigma %*% U
  UXsU <- t(U) %*% Xsigma %*% U
  UXsUinv <- solve(UXsU)
  sigma.X2 <- UtXsUt - UtXsU %*% UXsUinv %*% t(UtXsU)
  sigma.X2 <- (sigma.X2 + t(sigma.X2)) / 2


  X1 <- matrix(nrow = nrank, ncol = n, rnorm(n * nrank))
  mean.X2 <- UtXsU %*% UXsUinv %*% X1
  X2 <- mean.X2 + t(mvrnorm(ncol(mean.X2), rep(0, nrow(mean.X2)), sigma.X2))
  X <- t(solve(t(P)) %*% rbind(X1, X2))

  E <- matrix(nrow = n, ncol = q, rnorm(n * q, 0, 1))
  C <- U %*% D %*% t(V)
  Y3 <- X %*% U[, nrank] %*% t(V[, nrank]) * D[nrank, nrank]
  sigma <- sqrt(sum(Y3 ^ 2) / sum(E ^ 2)) / snr
  E <- E * sigma
  Y <- X %*% C + E

  return(list(Y = Y, X = X))
}