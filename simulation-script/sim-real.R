require(peer)
require(secure)
require(ggplot2)

setwd('d:/result')

data("CellCycle")
X <- CellCycle$X; Y <- CellCycle$Y 

peer.fit <- peer(X, Y, nrank = 4) 
d <- peer.fit$d; U <- peer.fit$U; V <- peer.fit$V
VD <- V %*% diag(d[1:ncol(V)])
factor.name <- c('First factor', 'Second factor', 'Third factor', 'Fourth factor') 

vd.table <- data.frame()
for (i in c(1:ncol(V))) {
  temp <- data.frame(loading = VD[, i], factor = rep(factor.name[i], nrow(V)), time = seq(0, 119, length = 18))
  vd.table <- rbind(vd.table, temp)
}

v.table <- data.frame()
for (i in c(1:ncol(V))) {
  temp <- data.frame(loading = V[, i], factor = rep(factor.name[i], nrow(V)), time = seq(0, 119, length = 18))
  v.table <- rbind(v.table, temp)
}

# print genes (peer)

genes <- colnames(X)

cat('genes (peer):\n')
for (i in c(1:ncol(U))) {
  ind <- U[, i] != 0 
  cat(factor.name[i], ": ", genes[ind], '\n')
}

# print genes (secure) 

secure.fit <- secure.path(Y, X, nrank = 4)
U <- secure.fit$U

cat('genes (secure):\n')
for (i in c(1:ncol(U))) {
  ind <- U[, i] != 0 
  cat(factor.name[i], ": ", genes[ind], '\n')
}
