library(scales)

#### Measurement Error

a.vals <- seq(6, 14, by = 0.04)
filenames <- list.files(path = "~/Dropbox/Projects/ERC-EPE/Output/sim_1", full.names = TRUE)
plotnames <- c("No Measurement Error",
               "Prediction Error but No Aggregation Error",
               "Aggregation Error but No Prediction Error",
               "Both Classical and Prediction Error")
idx <- c(1, 3, 7, 9)

pdf(file = "~/Dropbox/Projects/ERC-EPE/Output/plot_3.pdf", width = 9, height = 9)
par(mfrow = c(2,2))
j <- 1

for (k in idx){
  
  load(file = filenames[[k]])
  
  plot(a.vals, rslt$est[1,], type = "l", col = hue_pal()(6)[1], lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[j],
       ylim = c(0,0.1))
  grid(lty = 1)
  lines(a.vals, rslt$est[2,], type = "l", col = hue_pal()(6)[2], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[3,], type = "l", col = hue_pal()(6)[3], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[5,], type = "l", col = hue_pal()(6)[5], lwd = 2, lty = 1)
  
  if (k == idx[length(idx)]){
    
    legend(6, 0.1, legend=c("True ERC", "Naive", "RC", "BART+LOESS"),
           lty = hue_pal()(6)[c(1,2,3,5)], lwd=2, cex=0.8)
    
    
  }
  
  j <- j + 1
  
}

dev.off()

# Summary Table

tbl <- matrix(NA, nrow = length(rslt$est), ncol = 12)

for (k in 1:length(rslt$est)){
  
  bias <- round(colMeans(t(rslt$est[[k]]$bias)/rslt$est[[k]]$est[1,]), 3)
  se <- round(rowMeans(rslt$est[[k]]$se/rslt$est[[k]]$sd), 3)
  ci <- round(rowMeans(rslt$est[[k]]$cp), 3)
  
  tbl[k,] <- c(bias, se, ci)
  
}

colnames(tbl) <- outer(names(bias), c("Bias", "SE", "CI"), FUN = "paste")[1:12]
final <- cbind(rslt$scen_idx, tbl)

write.csv(final, file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/table_1.csv")