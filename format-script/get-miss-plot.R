require(ggplot2)
require(ggpubr)
setwd("d:/Git/PEER/Code/result")
path <- "./sim-miss"

snr <- 0.5


for (p in c(200,400)) {
  file = paste0(path, "/snr-", snr * 100, "-p-", p, ".RData")
  load(file)
  table$ErXC <- table$ErXC 
  table$miss <- as.character(table$miss)
  
  # ErC
  fig.c <- ggplot(table) +
    aes(x = miss, y = ErC) +
    geom_boxplot(notch = T, fill = "gray50") +
    xlab("missing rate") +
    ylab(expression("Er(C)" %*% 10 ^ 3)) +
    ylim(c(0,10))
  
  # ErXC
  fig.xc <- ggplot(table) +
    aes(x = miss, y = ErXC) +
    geom_boxplot(notch = T, fill = "gray50") +
    xlab("missing rate") +
    ylab("Er(XC)") +
    ylim(c(0,2))
  
  # FPR
  fig.fpr <- ggplot(table) +
    aes(x = miss, y = FPR) +
    geom_boxplot(notch = T, fill = "gray50") +
    xlab("missing rate") +
    ylab("FPR(%)") +
    ylim(c(0,11))
  
  # FNR
  fig.fnr <- ggplot(table) +
    aes(x = miss, y = FNR) +
    geom_boxplot(notch = T, fill = "gray50") +
    xlab("missing rate") +
    ylab("FNR(%)") +
    ylim(c(0,20))
  
  
  fig <-
    ggarrange(
      fig.c,
      fig.xc,
      fig.fpr,
      fig.fnr,
      ncol = 4,
      nrow = 1,
      common.legend = TRUE,
      legend = "bottom"
    )
  
  plot(fig)
  
  # file <- paste0("D:\\snr-", snr * 100, "-p-", p, ".eps")
  
  # ggsave(file, plot = fig, width = 7, height = 2, units = "cm")
  
  # setEPS()
  # postscript(paste0(file,".eps"),width = 2.75, height = 0.75)
  # plot(fig)
  # dev.off()
}
