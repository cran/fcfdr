## ---- include = FALSE,eval=FALSE----------------------------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>",
#    dev = 'png'
#  )
#  
#  Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  library(fcfdr)
#  library(cowplot)
#  library(ggplot2)
#  library(dplyr)
#  
#  data(T1D_application_data, package = "fcfdr")
#  head(T1D_application_data)

## ----eval=FALSE---------------------------------------------------------------
#  orig_p <- T1D_application_data$T1D_pval
#  chr <- T1D_application_data$CHR19
#  MAF <- T1D_application_data$MAF
#  q1 <- T1D_application_data$RA_pval
#  q2 <- T1D_application_data$DGF
#  q3 <- log(T1D_application_data$H3K27ac+1) # deal with long tail

## ----eval=FALSE---------------------------------------------------------------
#  ind_snps <- which(T1D_application_data$LDAK_weight != 0)

## ----eval=FALSE---------------------------------------------------------------
#  iter1_res <- flexible_cfdr(p = orig_p,
#                             q = q1,
#                             indep_index = ind_snps,
#                             maf = MAF)
#  
#  v1 <- iter1_res[[1]]$v

## ----eval=FALSE---------------------------------------------------------------
#  res1 <- data.frame(p = orig_p, q1, v1)
#  mid1 <- median(res1$q1)
#  
#  ggplot(res1, aes(x = p.adjust(p, method = "BH"), y = p.adjust(v1, method = "BH"), col = q1)) + geom_point(cex = 0.5) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_abline(intercept = 0, slope = 1,  linetype="dashed") + xlab("Original FDR") + ylab("V1 (FDR)") + ggtitle(paste0("Iteration 1")) + scale_color_gradient2(midpoint = mid1, low = "blue", mid = "white", high = "red", space = "Lab")

## ----eval=FALSE---------------------------------------------------------------
#  iter2_res <- binary_cfdr(p = v1,
#                           q = q2,
#                           group = chr)
#  
#  v2 <- iter2_res$v

## ----eval=FALSE---------------------------------------------------------------
#  res2 <- data.frame(p = v1, v2, q2)
#  res2$q2 <- as.factor(res2$q2)
#  
#  ggplot(res2, aes(x = p.adjust(p, method = "BH"), y = p.adjust(v2, method = "BH"), col = q2)) + geom_point(cex = 0.5) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_abline(intercept = 0, slope = 1,  linetype="dashed") + xlab("V1 (FDR)") + ylab("V2 (FDR)") + ggtitle(paste0("Iteration 2")) + scale_colour_manual(values = c("grey", "black"))

## ----eval=FALSE---------------------------------------------------------------
#  iter3_res <- flexible_cfdr(p = v2,
#                             q = q3,
#                             indep_index = ind_snps,
#                             maf = MAF)
#  v3 <- iter3_res[[1]]$v

## ----eval=FALSE---------------------------------------------------------------
#  res3 <- data.frame(p = v2, q3, v3)
#  
#  ggplot(res3, aes(x = p.adjust(p, method = "BH"), y = p.adjust(v3, method = "BH"), col = q3)) + geom_point(cex = 0.5) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_abline(intercept = 0, slope = 1,  linetype="dashed") + xlab("V2 (FDR)") + ylab("V3 (FDR)") + ggtitle(paste0("Iteration 3")) + scale_color_gradient2(midpoint = 1, low = "blue", mid = "white", high = "red", space = "Lab")

## ----eval=FALSE---------------------------------------------------------------
#  res <- data.frame(orig_p, q1 = iter1_res[[1]]$q, q2 = as.factor(iter2_res$q), q3 = iter3_res[[1]]$q, v1, v2, v3)
#  
#  head(res)

## ----eval=FALSE---------------------------------------------------------------
#  mid1 <- median(res$q1)
#  
#  ggplot(res, aes(x = p.adjust(orig_p, method = "BH"), y = p.adjust(v3, method = "BH"))) + geom_point(cex = 0.5, alpha = 0.5) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_abline(intercept = 0, slope = 1,  linetype="dashed", col = "red") + xlab("Original P (FDR)") + ylab("V3 (FDR)") + ggtitle(paste0("FDR adjusted v-values\nagainst original FDR values"))

## ----eval=FALSE---------------------------------------------------------------
#  ggplot(res, aes(x = -log10(orig_p), y = -log10(v3))) + geom_point(cex = 0.5, alpha = 0.5) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_abline(intercept = 0, slope = 1,  linetype="dashed", col = "red") + xlab("Original P (FDR) (-log10)") + ylab("V3 (FDR) (-log10)") + ggtitle(paste0("FDR adjusted v-values against\noriginal FDR values (FDR)")) + coord_cartesian(ylim = c(0,10), xlim = c(0,10))

## ----eval=FALSE---------------------------------------------------------------
#  p_fdr <- p.adjust(orig_p, method = "BH")
#  v3_fdr <- p.adjust(v3, method = "BH")
#  
#  # choose fdr threshold corresponding to genome-wide significance threshold
#  fdr_thr <- max(p_fdr[which(orig_p <= 5*10^{-8})])
#  
#  median(T1D_application_data$RA_p[which(v3_fdr < fdr_thr & p_fdr > fdr_thr)])
#  median(T1D_application_data$RA_p)
#  
#  mean(T1D_application_data$DGF[which(v3_fdr < fdr_thr & p_fdr > fdr_thr)])
#  mean(T1D_application_data$DGF)
#  
#  median(T1D_application_data$H3K27ac[which(v3_fdr < fdr_thr & p_fdr > fdr_thr)])
#  median(T1D_application_data$H3K27ac)

## ---- eval = FALSE------------------------------------------------------------
#  T1D_application_data$v3_fdr <- v3_fdr
#  
#  nCHR <- length(unique(T1D_application_data$CHR19))
#  T1D_application_data$BPcum <- NA
#  s <- 0
#  nbp <- c()
#  T1D_application_data <- data.frame(T1D_application_data)
#  for (i in unique(T1D_application_data$CHR19)){
#    nbp[i] <- max(T1D_application_data[T1D_application_data$CHR19 == i,]$BP19)
#    T1D_application_data[T1D_application_data$CHR19 == i,"BPcum"] <- T1D_application_data[T1D_application_data$CHR19 == i,"BP19"] + s
#    s <- s + nbp[i]
#  }
#  
#  axis.set <- T1D_application_data %>%
#    group_by(CHR19) %>%
#    summarize(center = (max(BPcum) + min(BPcum)) / 2)
#  
#  ggplot(T1D_application_data, aes(x = BPcum, y = -log10(v3_fdr), col = as.factor(CHR19))) + geom_point(cex = 0.75) + theme_cowplot(12) + background_grid(major = "xy", minor = "none") + geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed") + xlab("Position") +  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
#    scale_x_continuous(label = axis.set$CHR19, breaks = axis.set$center) + theme(legend.position = "none")+ theme(axis.text.x = element_text(size = 6, angle = 0)) + coord_cartesian(ylim=c(0,10)) + ylab(expression(paste("-log"[10],"(FDR)")))

