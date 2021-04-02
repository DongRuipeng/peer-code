require(peer)
require(secure)
require(doParallel)
require(foreach)

setwd("d:/result")
path <- "./sim-table-secure"
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

sim.setup.secure <- function(n = 100, p = 100, q = 100, nrank = 3, su = 4, sv = 5, snr = 0.5, miss = 0) {
  U <- matrix(0, p, nrank)
  V <- matrix(0, q, nrank)
  for (k in c(1:nrank)) {
    U[, k] <- c(rep(0, (k - 1) * su), sample(c(1, -1), su, replace = TRUE), rep(0, p - k * su))
    V[, k] <- c(rep(0, (k - 1) * sv), sample(c(1, -1), sv, replace = TRUE) * runif(sv, 0.3, 1), rep(0, q - k * sv))
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

sim.table.method2 <- function(control = control) {
  n <- control$n
  p <- control$p
  q <- control$q
  nrank <- control$nrank
  mrank <- control$mrank
  miss <- control$miss
  snr <- control$snr
  
  su <- 4
  sv <- 5
  sim.data <- sim.setup.secure(n, p, q, nrank, su, sv, snr, miss)
  X <- sim.data$X
  Y <- sim.data$Y
  U <- sim.data$U
  V <- sim.data$V
  d <- sim.data$d
  C <- sim.data$C
  
  secure.fit <- peer.fit <- data.frame(ErC = 0, ErXC = 0, FPR = 0, FNR = 0, time = 0, erank = 0)
  
  #secure(D is diagonal matrix decrease)
  start_time <- Sys.time()
  fit <- secure.path(Y, X, nrank = 4, nlambda = 100)
  end_time <- Sys.time()
  Chat <- fit$C.est
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- diag(fit$D)
  if (sum(dhat) == 0) {
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
  
  table <- rbind(secure.fit, peer.fit)
  method <- c("SeCURE", "PEER")
  table <- cbind(table, method = method)
  rownames(table) <- NULL
  return(table)
}

n.fix <- 200
q.fix <- 100
r.fix <- 3
miss.fix <- 0.1

cl.num <- 1
repetition <- 200

sink(paste0(path, "/log"))

cl <- makeCluster(cl.num)
registerDoParallel(cl)

for (snr in c(0.5,1)) {
  cat('begin snr =', snr, '\n')
  t <- Sys.time()
  for (p in c(200,400)) {
    # define simulation setup
    setting <- sim.control(n = n.fix, p = p, q = q.fix, nrank = r.fix, snr = snr, miss = miss.fix)
    # begin simulation
    table <- foreach(i = 1:repetition, .combine = 'rbind', .packages = c('peer', 'secure'), .errorhandling = 'pass') %dopar% {
      sim.table.method2(control = setting)
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
