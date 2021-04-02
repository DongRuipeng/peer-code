setwd("d:/Git/PEER/Code/result")
require(xlsx)

my.mean <- function(x) {
  return(mean(x, na.rm = TRUE))
}

my.sd <- function(x) {
  return(sd(x, na.rm = TRUE))
}

path <- "./sim-table-v-3"
# sink(paste0(path,"/table-default.tex"))

miss.fix <- 0.1

for (snr in c(0.25, 0.5)) {
  for (p in c(200, 400)) {
    load(paste0(path, "/miss-", miss.fix * 100, "-snr-", snr * 100, "-p-", p, ".Rdata"))

    table$time <- as.double(table$time)
    performance <- colnames(table)[1:6]
    ave <- std <- data.frame()
    for (method in c("mRRR","SeCURE", "PEER")) {
      temp <- table[table$method == method, 1:length(performance)]
      colnames(temp) <- NULL
      ave <- rbind(ave, apply(temp, 2, my.mean))
      std <- rbind(std, apply(temp, 2, my.sd))
    }
    method.name <- c("mRRR","SeCURE", "PEER")
    rownames(ave) <- method.name
    rownames(std) <- method.name
    colnames(ave) <- performance
    colnames(std) <- performance

    # cat('p =', p, 'snr =', snr, '\n\n')
    # cat('Ave: \n')
    # print(ave)
    # cat('\n')
    # cat('Std: \n')
    # print(std)
    # cat('\n\n')
    
    filename <- paste0(path, "/snr-", snr * 100, "-p-", p, ".xlsx")
    
    write.xlsx(ave, file = filename, sheetName = "ave")
    write.xlsx(std,file = filename,sheetName = "std",append = TRUE)
  }
}

# sink()