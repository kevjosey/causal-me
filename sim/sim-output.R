library(ggplot2)
library(cowplot)
library(ggpubr)

#### Measurement Error

a.vals <- seq(6, 14, by = 0.04)
filenames_1 <- c("~/Dropbox/Projects/ERC-EPE/Output/sim_1/800_5_0_0_a_a_a.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_1/800_5_0_1_a_a_a.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_1/800_5_1.4142135623731_0_a_a_a.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_1/800_5_1.4142135623731_1_a_a_a.RData")
plotnames <- c("No Measurement Error",
               "Prediction Error/No Aggregation Error",
               "Aggregation Error/No Prediction Error",
               "Aggregation Error and Prediction Error")
# idx <- c(1,2,7,8)
# idx <- c(10,11,16,17)
# idx <- c(19,20,25,26)
# idx <- c(39,40,45,46)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_list <- list()
j <- 1

for (fdx in filenames_1){
  
  load(file = fdx)
  
  plot_dat <- data.frame(estimate = c(unname(rslt$est[1,]), unname(rslt$est[2,]) ,unname(rslt$est[3,]),unname(rslt$est[4,]),unname(rslt$est[5,])),
                         label = factor(rep(c("True ERF", "No Correction", "Regression Calibration","BART Multiple Imputation", "GLM Multiple Imputation"), each = 201),
                                        levels = c("True ERF", "No Correction", "Regression Calibration", "BART Multiple Imputation", "GLM Multiple Imputation")),
                         exposure = as.numeric(colnames(rslt$est, 5)))
  
  mse_dat <- data.frame(mse = rowMeans(sqrt(rslt$mse)), 
                        x = c(7, 9, 11, 13),
                        label = factor(c("No Correction", "Regression Calibration", "BART Multiple Imputation", "GLM Multiple Imputation"),
                                       levels = c("No Correction", "Regression Calibration", "BART Multiple Imputation", "GLM Multiple Imputation")))
  
  erf_plot <- plot_dat %>% 
    ggplot(aes(x = exposure, y = estimate, color = label)) + 
    geom_line(size = 1) + 
    scale_color_manual(values=c(gg_color_hue(5))) +
    geom_bar(mse_dat, mapping = aes(fill = label, x = x, y = mse*10), stat = "identity") +
    scale_fill_manual(values = alpha(c(gg_color_hue(5)[2:5]), 0.3))+
    labs(x = "Exposure", y = "Response Rate", color = "", fill = "") + 
    scale_y_continuous(limits = c(0,0.3), sec.axis = sec_axis(~ ./10 , name = "Average Root Mean Squared Error"))+
    guides(fill = "none") +
    theme_bw() + ggtitle(plotnames[j]) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="none")

  
  plot_list[[j]] <- erf_plot
  j <- j + 1
  
}

erf_plot <- plot_dat %>% 
  ggplot(aes(x = exposure, y = estimate, color = label)) + 
  geom_line(size = 1) + 
  scale_color_manual(values=c(gg_color_hue(5))) +
  labs(x = "Exposure", y = "Response Rate", color = "", fill = "") + 
  theme_bw() + ggtitle(plotnames[j]) +
  theme(plot.title = element_text(hjust = 0.5))

legend_b <- get_legend(erf_plot + theme(legend.position="bottom"))

pdf(file = "~/Dropbox/Projects/ERC-EPE/Output/plot_sim.pdf", width = 9, height = 9, onefile = FALSE)
ggarrange(plotlist = plot_list, ncol = 2, nrow = 2, legend = "bottom", common.legend = TRUE, legend.grob = legend_b)
dev.off()

### Model Misspecification

a.vals <- seq(6, 14, by = 0.04)

filenames_2 <- c("~/Dropbox/Projects/ERC-EPE/Output/sim_2/800_5_1.4142135623731_1_b_a_a.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_2/800_5_1.4142135623731_1_a_b_a.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_2/800_5_1.4142135623731_1_a_a_b.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_2/800_5_1.4142135623731_1_a_b_b.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_2/800_5_1.4142135623731_1_b_a_b.RData",
                 "~/Dropbox/Projects/ERC-EPE/Output/sim_2/800_5_1.4142135623731_1_b_b_a.RData")

plotnames <- c("Incorrect GPS/Correct Outcome/Correct EPE",
               "Correct GPS/Incorrect Outcome/Correct EPE",
               "Correct GPS/Correct Outcome/Incorrect EPE",
               "Correct GPS/Incorrect Outcome/Incorrect EPE",
               "Incorrect GPS/Correct Outcome/Incorrect EPE",
               "Incorrect GPS/Incorrect Outcome/Correct EPE")

pdf(file = "~/Dropbox/Projects/ERC-EPE/Output/plot_mis.pdf", width = 8, height = 12)
par(mfrow = c(3,2))
j <- 1

for (fdx in filenames_2){
  
  load(file = fdx)
  
  plot(a.vals, rslt$est[1,], type = "l", col = "black", lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[j],
       ylim = c(0,0.3))
  grid(lty = 1)
  lines(a.vals, rslt$est[5,], type = "l", col = gg_color_hue(5)[4], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[6,], type = "l", col = gg_color_hue(5)[5], lwd = 2, lty = 1)
  lines(a.vals, rslt$upper[4,], type = "l", col = gg_color_hue(5)[4], lwd = 2, lty = 2)
  lines(a.vals, rslt$upper[5,], type = "l", col = gg_color_hue(5)[5], lwd = 2, lty = 2)
  lines(a.vals, rslt$lower[4,], type = "l", col = gg_color_hue(5)[4], lwd = 2, lty = 2)
  lines(a.vals, rslt$lower[5,], type = "l", col = gg_color_hue(5)[5], lwd = 2, lty = 2)
  
  if (j == 6){
    
    legend(6, 0.3, legend=c("True ERF", "GLM", "DR", "95% CI"),
           col = c("black", gg_color_hue(4)[4:3], "black"),
           bg = "white", lty = c(1,1,1,2), lwd=2, cex=0.8)
    
    
  }
  
  j <- j + 1
  
}

dev.off()

### Summary Table

tbl <- data.frame()
filenames <- c(list.files(path = "~/Dropbox/Projects/ERC-EPE/Output/sim_1", full.names = TRUE),
               list.files(path = "~/Dropbox/Projects/ERC-EPE/Output/sim_2", full.names = TRUE))
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
write.csv(tbl, file = "~/Dropbox/Projects/ERC-EPE/Output/table.csv")
