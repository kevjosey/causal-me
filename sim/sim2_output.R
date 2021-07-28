
library(scales)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.04)
filenames <- list.files(path = "~/Dropbox/Projects/ERC-EPE/Output/sim2", full.names = TRUE)

# Summary Plot

plotnames <- c("GPS: \"a\"; Outcome: \"a\"",
               "GPS: \"b\"; Outcome: \"a\"",
               "GPS: \"a\"; Outcome: \"b\"",
               "GPS: \"b\"; Outcome: \"b\"")
idx <- c(22,23,24)

pdf(file = "~/Dropbox/Projects/ERC-EPE/Output/plot_2.pdf", width = 9, height = 9)
par(mfrow = c(2,2))

for (k in idx){
  
  load(file = filenames[[k]])
  
  plot(a.vals, rslt$est[1,], type = "l", col = hue_pal()(6)[1], lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[k],
       ylim = c(0,0.15))
  grid(lty = 1)
  lines(a.vals, rslt$est[2,], type = "l", col = hue_pal()(6)[2], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[3,], type = "l", col = hue_pal()(6)[3], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[4,], type = "l", col = hue_pal()(6)[4], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[5,], type = "l", col = hue_pal()(6)[5], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[6,], type = "l", col = hue_pal()(6)[6], lwd = 2, lty = 1)
  
  if (k == idx[length(idx)]){
    
    legend(6, 0.15, legend=c("True ERF", "Observed", "Naive", "BLP", "BART"),
           col=hue_pal()(5),
           lty = c(1,1,1,1,1), lwd=2, cex=0.8)
    
    
  }
  
}

dev.off()

# Summary Table

tbl <- data.frame()

for (k in 1:length(filenames)) {
  
  load(filenames[k])
  
  tbl[k,1:6] <- rslt$sceario
  
  bias <- round(colMeans(t(rslt$bias)/rslt$est[1,]), 3)
  se <- round(rowMeans(rslt$se/rslt$sd), 3)
  ci <- round(rowMeans(rslt$cp), 3)
  
  tbl[k,7:21] <- c(bias, se, ci)
  
  
}

colnames(tbl)[1:6] <- colnames(rslt$sceario)
colnames(tbl)[7:21] <- outer(rownames(rslt$bias), c("Bias", "SE", "CI"), FUN = "paste")[1:15]

write.csv(tbl, file = "~/Dropbox/Projects/ERC-EPE/Output/table_2.csv")