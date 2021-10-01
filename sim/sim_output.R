#### Measurement Error

a.vals <- seq(6, 14, by = 0.04)
filenames_1 <- list.files(path = "~/Dropbox/Projects/ERF-EPE/Output/sim_1", full.names = TRUE)
plotnames <- c("No Measurement Error",
               "Prediction Error but No Aggregation Error",
               "Aggregation Error but No Prediction Error",
               "Both Aggregation and Prediction Error")
# idx <- c(1,2,7,8)
# idx <- c(10,11,16,17)
# idx <- c(19,20,25,26)
idx <- c(28,29,34,35)

pdf(file = "~/Dropbox/Projects/ERF-EPE/Output/plot_sim.pdf", width = 9, height = 9)
par(mfrow = c(2,2))
j <- 1

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

for (k in idx){
  
  load(file = filenames_1[[k]])
  
  plot(a.vals, rslt$est[1,], type = "l", col = "black", lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[j],
       ylim = c(0,0.1))
  grid(lty = 1)
  lines(a.vals, rslt$est[2,], type = "l", col = gg_color_hue(3)[1], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[3,], type = "l", col = gg_color_hue(3)[2], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[5,], type = "l", col = gg_color_hue(3)[3], lwd = 2, lty = 1)
  
  if (k == idx[length(idx)]){
    
    legend(6, 0.1, legend=c("True ERF", "No Correction", "Regression Calibration", "Multiple Imputation"),
           col = c("black", gg_color_hue(3)), lwd=2, cex=0.8)
    
    
  }
  
  j <- j + 1
  
}

dev.off()

### Model Misspecification

a.vals <- seq(6, 14, by = 0.04)
filenames_2 <- list.files(path = "~/Dropbox/Projects/ERF-EPE/Output/sim_2", full.names = TRUE)

plotnames <- c("GPS: \"a\"; Outcome: \"a\"; EPE: \"b\"",
               "GPS: \"a\"; Outcome: \"b\"; EPE: \"a\"",
               "GPS: \"a\"; Outcome: \"b\"; EPE: \"b\"",
               "GPS: \"b\"; Outcome: \"a\"; EPE: \"a\"",
               "GPS: \"b\"; Outcome: \"a\"; EPE: \"b\"",
               "GPS: \"b\"; Outcome: \"b\"; EPE: \"a\"")
idx <- c(2:7)

pdf(file = "~/Dropbox/Projects/ERF-EPE/Output/plot_mis.pdf", width = 11, height = 11)
par(mfrow = c(3,2))
j <- 1

for (k in idx){
  
  load(file = filenames_2[[k]])
  
  plot(a.vals, rslt$est[1,], type = "l", col = "black", lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[j],
       ylim = c(0,0.1))
  grid(lty = 1)
  lines(a.vals, rslt$est[6,], type = "l", col = gg_color_hue(4)[4], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[5,], type = "l", col = gg_color_hue(4)[3], lwd = 2, lty = 1)
  lines(a.vals, rslt$upper[5,], type = "l", col = gg_color_hue(4)[4], lwd = 2, lty = 2)
  lines(a.vals, rslt$upper[4,], type = "l", col = gg_color_hue(4)[3], lwd = 2, lty = 2)
  lines(a.vals, rslt$lower[5,], type = "l", col = gg_color_hue(4)[4], lwd = 2, lty = 2)
  lines(a.vals, rslt$lower[4,], type = "l", col = gg_color_hue(4)[3], lwd = 2, lty = 2)
  
  if (k == idx[length(idx)]){
    
    legend(6, 0.1, legend=c("True ERF", "GLM", "BART", "95% CI"),
           col = c("black", gg_color_hue(4)[4:3], "black"),
           lty = c(1,1,1,2), lwd=2, cex=0.8)
    
    
  }
  
  j <- j + 1
  
}

dev.off()

### Summary Table

tbl <- data.frame()
filenames <- c(unlist(filenames_1), unlist(filenames_2))

for (k in 1:length(filenames)){
  
  load(file = filenames[[k]])
  
  idx <- 126
  
  relvec <- rslt$est[1,]
  relmat <- matrix(rep(relvec, times = 4), nrow = 4, byrow = TRUE)
  
  tbl[k,1:7] <- rslt$scenario
  tbl[k,8:11] <- paste(round(colMeans(t(rslt$est[2:5,] - relmat)/relvec), 2), " (", round(c(rslt$est[2:5,idx] - relmat[,idx])/relvec[idx], 2), ")", sep = "")
  tbl[k,12:15] <- paste(round(rowMeans(sqrt(rslt$mse[1:4,])), 3), " (", round(sqrt(rslt$mse[1:4,idx]), 3), ")", sep = "")
  tbl[k,16:19] <- paste(round(rowMeans(rslt$cp[1:4,]), 2), " (", round(rslt$cp[1:4,idx], 2), ")", sep = "")
  
}

colnames(tbl)[8:19] <- outer(rownames(rslt$bias[1:4,]), c("Bias", "MSE", "CI"), FUN = "paste")[1:12]

write.csv(tbl, file = "~/Dropbox/Projects/ERF-EPE/Output/table.csv")
