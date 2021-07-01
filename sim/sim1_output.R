# Summary Plot

load(file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/sim_1.RData")
plotnames <- c("No Measurement Error",
               "Classical Error but No Prediction Error",
               "Prediction Error but No Classical Error",
               "Both Classical and Prediction Error")
idx <- c(20,24,28,32)

filename <- paste0("~/Dropbox (Personal)/Projects/ERC-EPE/Output/plot_1.pdf")
pdf(file = filename, width = 9, height = 9)
par(mfrow = c(2,2))

for (k in 1:4){
  
  plot(a.vals, rslt$est[[idx[k]]]$est[1,], type = "l", col = "darkgreen", lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[k],
       ylim = c(0,0.08))
  grid(lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[2,], type = "l", col = "red", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[idx[k]]]$est[3,], type = "l", col = "blue", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[idx[k]]]$est[4,], type = "l", col = "red", lwd = 2, lty = 3)
  lines(a.vals, rslt$est[[idx[k]]]$est[5,], type = "l", col = "blue", lwd = 2, lty = 3)
  
  if (k == 4){
    
    legend(x = 7.7, y = 0.025, legend=c("True ERF", "Without Prediction Correction",
                                        "With Prediction Correction", "Without Classical Correction",
                                        "With Classical Correction"),
           col=c("darkgreen", "red", "blue", "black", "black"),
           lty = c(1,1,1,2,3), lwd=2, cex=0.8)
    
  }
  
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