setwd("/Users/dr/Desktop/MR.SMK-T2D/updated_forestplots")
pdf("cigarettes_per_day.CAD.forestplot.pdf")
res_single <- mr_singlesnp(cad_cig_day_dat, all_method = c("mr_ivw",
                                               "mr_egger_regression",
                                               "mr_weighted_median"))
singlesnp_results <- res_single
exponentiate <- FALSE

requireNamespace("ggplot2", quietly = TRUE)
requireNamespace("plyr", quietly = TRUE)
res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d) {
  d <- plyr::mutate(d)
  if (sum(!grepl("All", d$SNP)) < 2) {
    return(blank_plot("Insufficient number of SNPs"))
  }
  levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "Inverse variance weighted"
  levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "MR Egger"
  levels(d$SNP)[levels(d$SNP) == "All - Weighted median"] <- "Weighted median"
  #d[d$SNP == "All - Inverse variance weighted", "SNP"] <- "Inverse variance weighted"
  #d[d$SNP == "All - MR Egger", "SNP"] <- "MR Egger"
  #d[d$SNP == "All - Weighted median", "SNP"] <- "Inverse variance weighted"
                     am <- grep("All", d$SNP, value = TRUE)

                     d$up <- d$b + 1.96 * d$se
                     d$lo <- d$b - 1.96 * d$se
                     
                     # change unit
                     # binary: initiation + cessation
                     #d$b <- exp(d$b * log(2))
                     #d$up <- exp(d$up * log(2))
                     #d$lo <- exp(d$lo * log(2))
                     
                     # continuous: cig_day
                     d$b <- exp(d$b)
                     d$up <- exp(d$up)
                     d$lo <- exp(d$lo)
                     
                     d$tot <- 0.01
                     d$tot[d$SNP %in% am] <- 1
                     d$SNP <- as.character(d$SNP)
                     nom <- d$SNP[!d$SNP %in% am]
                     nom <- nom[order(d$b)]
                     d <- rbind(d, d[nrow(d), ])
                     #d$SNP[nrow(d) - 1] <- ""
                     #d$b[nrow(d) - 1] <- NA
                     #d$up[nrow(d) - 1] <- NA
                     #d$lo[nrow(d) - 1] <- NA
                     d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
                     xint <- 0
                     if (exponentiate) {
                       d$b <- exp(d$b)
                       d$up <- exp(d$up)
                       d$lo <- exp(d$lo)
                       xint <- 1
                     }
                     #print(tail(d, 4))
                     d <- tail(d, 4)
                     d <- head(d, 3)
                     d[d$SNP == "Inverse variance weighted", "samplesize"] <- 3
                     d[d$SNP == "Weighted median", "samplesize"] <- 2
                     d[d$SNP == "MR Egger", "samplesize"] <- 1
                     d$SNP2 <- reorder(d$SNP, d$samplesize)
                     print(d)
                     ggplot2::ggplot(d, aes(y = SNP2, x = b)) + 
                       ggplot2::geom_vline(xintercept = 1, linetype = "dotted") + 
                       ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, size=as.factor(tot),colour = as.factor(tot)), size=4,
                                height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot)), size=16)  + 
                       ggplot2::scale_colour_manual(values = c("black", "red")) + 
                       ggplot2::scale_size_manual(values = c(0.3, 1)) + 
                       ggplot2::theme(legend.position = "none", 
                        axis.text.y = ggplot2::element_text(size = 14), 
                        axis.ticks.y = ggplot2::element_line(size = 0), 
                        axis.title.x = ggplot2::element_text(size = 14)) + 
                       ggplot2::labs(y = "", x = paste0("MR odds ratio for\n'", 
                                                        d$exposure[1], "' on '", d$outcome[1], "'"))
                   })
p2 <- res

p2[[1]]
dev.off()