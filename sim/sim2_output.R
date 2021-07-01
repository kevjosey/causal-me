# Summary Plot

load(file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/sim_2.RData")
plotnames <- c("GPS: \"a\"; Outcome: \"a\"",
               "GPS: \"b\"; Outcome: \"a\"",
               "GPS: \"a\"; Outcome: \"b\"",
               "GPS: \"b\"; Outcome: \"b\"")
idx <- c(4,8,12,16)

filename <- paste0("~/Dropbox (Personal)/Projects/ERC-EPE/Output/plot_2.pdf")
pdf(file = filename, width = 9, height = 9)
par(mfrow = c(2,2))

for (k in 1:4){
  
  plot(a.vals, rslt$est[[idx[k]]]$est[1,], type = "l", col = hue_pal()(6)[1], lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[k],
       ylim = c(0,0.12))
  grid(lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[2,], type = "l", col = hue_pal()(6)[2], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[3,], type = "l", col = hue_pal()(6)[3], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[4,], type = "l", col = hue_pal()(6)[4], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[5,], type = "l", col = hue_pal()(6)[5], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[5,], type = "l", col = hue_pal()(6)[6], lwd = 2, lty = 1)
  
  if (k == 4){
    
    legend(6, 0.15, legend=c("True ERF", "Observed", "Naive", "BLP", "BART"),
           col=hue_pal()(5),
           lty = c(1,1,1,1,1), lwd=2, cex=0.8)
    
    
  }
  
}

dev.off()

# Summary Table

tbl <- matrix(NA, nrow = length(rslt$est), ncol = 6)

for (k in 1:length(rslt$est)){
  
  bias <- round(colMeans(t(rslt$est[[k]]$bias)/rslt$est[[k]]$est[1,]), 3)
  se <- round(rowMeans(rslt$est[[k]]$se/rslt$est[[k]]$sd), 3)
  ci <- round(rowMeans(rslt$est[[k]]$cp), 3)
  
  tbl[k,] <- c(bias, se, ci)
  
}

colnames(tbl) <- outer(names(bias), c("Bias", "SE", "CI"), FUN = "paste")[1:6]
final <- cbind(rslt$scen_idx, tbl)

write.csv(final, file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/table_2.csv")