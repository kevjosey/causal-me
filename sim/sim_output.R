library(scales)

#### Measurement Error

a.vals <- seq(6, 14, by = 0.04)
filenames <- list.files(path = "~/Dropbox/Projects/ERC-EPE/Output/sim_1", full.names = TRUE)
plotnames <- c("No Measurement Error",
               "Prediction Error but No Aggregation Error",
               "Aggregation Error but No Prediction Error",
               "Both Classical and Prediction Error")
idx <- c(1,2,7,8)

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
    
    legend(6, 0.1, legend=c("True ERF", "Naive", "Regression Calibration", "Multiple Imputation"),
           col = hue_pal()(6)[c(1,2,3,5)], lwd=2, cex=0.8)
    
    
  }
  
  j <- j + 1
  
}

dev.off()

# Summary Table

tbl <- data.frame()

for (k in 1:length(filenames)){
  
  load(file = filenames[[k]])
  
  tbl[k,1:7] <- rslt$scenario
  tbl[k,8:11] <- round(colMeans(t(rslt$bias)/rslt$est[1,]), 3)
  tbl[k,12:15] <- round(rowMeans(rslt$se), 3)
  tbl[k,16:19] <- round(rowMeans(rslt$sd), 3)
  tbl[k,20:23] <- round(rowMeans(rslt$cp), 3)
  
}

colnames(tbl)[8:23] <- outer(names(bias), c("Bias", "SE", "SD", "CI"), FUN = "paste")[1:16]
final <- cbind(rslt$scen_idx, tbl)

write.csv(final, file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/table_1.csv")

# Summary Plot

plotnames <- c("GPS: \"A\"; Outcome: \"A\"; EPE: \"B\"",
               "GPS: \"A\"; Outcome: \"B\"; EPE: \"A\"",
               "GPS: \"A\"; Outcome: \"B\"; EPE: \"B\"",
               "GPS: \"B\"; Outcome: \"A\"; EPE: \"A\"",
               "GPS: \"B\"; Outcome: \"A\"; EPE: \"B\"",
               "GPS: \"B\"; Outcome: \"B\"; EPE: \"A\"")
idx <- c(50:55)

pdf(file = "~/Dropbox/Projects/ERC-EPE/Output/plot_2.pdf", width = 11, height = 9)
par(mfrow = c(3,2))
j <- 1

for (k in idx){
  
  load(file = filenames[[k]])
  
  plot(a.vals, rslt$est[1,], type = "l", col = hue_pal()(3)[1], lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[j],
       ylim = c(0,0.15))
  grid(lty = 1)
  # lines(a.vals, rslt$est[2,], type = "l", col = hue_pal()(5)[2], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[3,], type = "l", col = hue_pal()(3)[2], lwd = 2, lty = 1)
  # lines(a.vals, rslt$est[4,], type = "l", col = hue_pal()(5)[4], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[5,], type = "l", col = hue_pal()(3)[3], lwd = 2, lty = 1)
  # lines(a.vals, rslt$est[2,] - 1.96*rslt$se[1,], type = "l", col = hue_pal()(5)[2], lwd = 2, lty = 2)
  lines(a.vals, rslt$est[3,] - 1.96*rslt$se[2,], type = "l", col = hue_pal()(3)[2], lwd = 2, lty = 2)
  # lines(a.vals, rslt$est[4,] - 1.96*rslt$se[3,], type = "l", col = hue_pal()(5)[4], lwd = 2, lty = 2)
  lines(a.vals, rslt$est[5,] - 1.96*rslt$se[4,], type = "l", col = hue_pal()(3)[3], lwd = 2, lty = 2)
  # lines(a.vals, rslt$est[2,] + 1.96*rslt$se[1,], type = "l", col = hue_pal()(5)[2], lwd = 2, lty = 2)
  lines(a.vals, rslt$est[3,] + 1.96*rslt$se[2,], type = "l", col = hue_pal()(3)[2], lwd = 2, lty = 2)
  # lines(a.vals, rslt$est[4,] + 1.96*rslt$se[3,], type = "l", col = hue_pal()(5)[4], lwd = 2, lty = 2)
  lines(a.vals, rslt$est[5,] + 1.96*rslt$se[4,], type = "l", col = hue_pal()(3)[3], lwd = 2, lty = 2)
  
  if (k == idx[length(idx)]){
    
    legend(6, 0.15, legend=c("True ERC", "RC", "BART+LOESS", "95% CI"),
           col= c(hue_pal()(3), "#000000"),
           lty = c(1,1,1,2), lwd=2, cex=0.8)
    
    
  }
  
  j <- j + 1
  
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