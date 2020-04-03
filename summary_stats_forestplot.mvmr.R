library(ggplot2)

# this is from mvmr final result
## %%% need to change df
m1 <- mvmr_out$coef
m1
data <- data.frame("beta" = m1[,"Estimate"], "se" = m1[,"Std. Error"], 
           df = 283, "t" = m1[,"t value"], "exposure" = c("Smoking", "BMI"))
data$up <- data$beta + data$se * qt(0.975, df = data$df)
data$low <- data$beta - data$se * qt(0.975, df = data$df)

# exp
data$exp.beta <- exp(data$beta)
data$exp.up <- exp(data$up)
data$exp.low <- exp(data$low)

data$exp2.beta <- exp(log(2) * data$beta)
data$exp2.up <- exp(log(2) * data$up)
data$exp2.low <- exp(log(2) * data$low)

# plot
## %%% need to change filename and title
## %%% need to change scale of beta
pdf("/Users/dr/Desktop/MR.SMK-T2D/better_instrument_forestplots/mvmr.smk-bmi->HbA1c.pdf")
ggplot2::ggplot(data, aes(y = exposure, x = exp2.beta)) + 
  ggplot2::geom_vline(xintercept = 1, linetype = "dotted") + 
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = exp2.low, xmax = exp2.up, size=as.factor(0.01),
                                       colour = as.factor(0.01)), size=4,
                          height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(0.01)), size=16)  + 
  ggplot2::scale_colour_manual(values = c("black", "red")) + 
  ggplot2::scale_size_manual(values = c(0.3, 1)) + 
  ggplot2::theme(legend.position = "none", 
                 axis.title.y = ggplot2::element_text(size = 14), 
                 axis.ticks.y = ggplot2::element_line(size = 0), 
                 axis.title.x = ggplot2::element_text(size = 14)) + 
  ggplot2::ggtitle("MR odds ratio for Smoking on HbA1c conditioning on BMI") + 
  theme(plot.title = element_text(size = 15, face = "bold"))
dev.off()



