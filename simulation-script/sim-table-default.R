require(peer)
require(rrpack)
require(secure)
require(doParallel)
require(foreach)

setwd("d:/result")
path <- "./sim-table-v-3"
dir.create(path = path, showWarnings = FALSE, recursive = TRUE)

peer.threshold <- function(X, Y, mrank = 20) {
  t_start <- Sys.time()
  erank <- est.rank(Y, mrank)
  t_end <- Sys.time()
  t <- difftime(t_end, t_start, units = "secs")
  fit <- peer(X, Y, nrank = erank, penalty = "l1")
  fit$time <- as.double(fit$time) + as.double(t)
  fit <- c(fit, erank = erank)
  return(fit)
}

cv.mrrr.rank <- function(X, Y, rank.list = c(2:6), nfold = 5) {
  n <- nrow(X)
  K <- floor(n / nfold)
  q <- ncol(Y)
  i <- 1
  error <- rep(0, length(rank.list))
  family <- list(gaussian())
  for (rank in rank.list) {
    error.temp <- rep(0, nfold)
    for (t in c(1:nfold)) {
      ind.test <- c((1 + K * (t - 1)):(K * t))
      X.train <- X[-ind.test,]
      Y.train <- Y[-ind.test,]
      X.test <- X[ind.test,]
      Y.test <- Y[ind.test,]

      fit <- mrrr(Y.train, X.train, family = family, penstr = list(penaltySVD = "rankCon", lambdaSVD = rank))
      Chat <- coef(fit)[-1,]
      A <- Y.test - X.test %*% Chat
      ind.na <- is.na(A)
      A[ind.na] <- 0
      error.temp[t] <- norm(A, 'F') / sqrt(K * q - sum(ind.na))
    }
    error[i] <- mean(error.temp)
    i <- i + 1
  }
  i_hat <- which.min(error)
  return(rank.list[i_hat])
}

sim.table.default <- function(control = control) {
  n <- control$n
  p <- control$p
  q <- control$q
  nrank <- control$nrank
  mrank <- control$mrank
  s <- control$s
  miss <- control$miss
  snr <- control$snr
  
  sim.data <- sim.setup(n, p, q, nrank, s, snr, miss)
  X <- sim.data$X
  Y <- sim.data$Y
  U <- sim.data$U
  V <- sim.data$V
  d <- sim.data$d
  C <- sim.data$C
  
  mrrr.fit <- secure.fit <- peer.fit <- data.frame(ErC = 0, ErXC = 0, FPR = 0, FNR = 0, time = 0, erank = 0)
  
  # mrrr 
  family <- list(gaussian())
  start_time <- Sys.time()
  erank <- cv.mrrr.rank(X, Y)
  fit <- mrrr(Y, X, family = family, penstr = list(penaltySVD = "rankCon", lambdaSVD = erank))
  end_time <- Sys.time()
  Chat <- coef(fit)[-1,]
  svdXC <- svd(X %*% Chat / sqrt(n), nu = erank, nv = erank)
  Vhat <- svdXC$v
  dUhat <- Chat %*% Vhat
  dhat <- svdXC$d[1:erank]
  Uhat <- dUhat %*% diag(dhat ^ (-1))
  mrrr.fit$erank <- erank 
  mrrr.fit$ErC <- est.err(C, Chat)
  mrrr.fit$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  mrrr.fit$FPR <- rate$fprate
  mrrr.fit$FNR <- rate$fnrate
  mrrr.fit$time <- difftime(end_time, start_time, units = "secs")
  
  #secure(D is diagonal matrix decrease)
  start_time <- Sys.time()
  fit <- secure.path(Y, X, nrank = 4, nlambda = 100)
  end_time <- Sys.time()
  Chat <- fit$C.est
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- diag(fit$D)
  if(sum(dhat) == 0) {
    erank <- 0
  } else {
    erank <- length(dhat)
  }
  secure.fit$erank <- erank 
  secure.fit$ErC <- est.err(C, Chat)
  secure.fit$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  secure.fit$FPR <- rate$fprate
  secure.fit$FNR <- rate$fnrate
  secure.fit$time <- difftime(end_time, start_time, units = "secs")
  
  # peer (d is a vector)
  fit <- peer.threshold(X, Y)
  Chat <- fit$C
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- fit$d
  peer.fit$erank <- fit$erank 
  peer.fit$ErC <- est.err(C, Chat)
  peer.fit$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  peer.fit$FPR <- rate$fprate
  peer.fit$FNR <- rate$fnrate
  peer.fit$time <- fit$time
  
  table <- rbind(mrrr.fit, secure.fit, peer.fit)
  method <- c("mRRR", "SeCURE", "PEER")
  table <- cbind(table, method = method)
  rownames(table) <- NULL
  return(table)
}

n.fix <- 100
q.fix <- 100
r.fix <- 3
s.fix <- 4
miss.fix <- 0.1

cl.num <- 1
repetition <- 200

sink(paste0(path, "/log"))

cl <- makeCluster(cl.num)
registerDoParallel(cl)

for (snr in c(0.25,0.5)) {
  cat('begin snr =', snr, '\n')
  t <- Sys.time()
  for (p in c(200,400)) {
    # define simulation setup
    setting <- sim.control(n = n.fix, p = p, q = q.fix, nrank = r.fix, s = s.fix, snr = snr, miss = miss.fix)
    # begin simulation
    table <- foreach(i = 1:repetition, .combine = 'rbind', .packages = c('peer', 'rrpack', 'secure'), .errorhandling = 'pass') %dopar% {
      sim.table.default(control = setting)
    }
    filename <- paste0(path, "/miss-", miss.fix * 100, "-snr-", snr * 100, "-p-", p, ".RData")
    save(table, file = filename)
  }
  cat('snr =', snr, 'is finished, elapsed time is', format(Sys.time() - t), '\n')
}

# end parallel computing
stopImplicitCluster()
stopCluster(cl)

sink()
