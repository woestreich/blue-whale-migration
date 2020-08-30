library(tidyverse)

rm(list = ls())

ci_daily <- read.csv("ci_daily.csv",header = FALSE)
colnames(ci_daily) <- c("CI","DCI","month")

ci_sig <- data.frame(as.numeric(rep(1,12)))
colnames(ci_sig) <- c("pv")
dci_sig <- data.frame(as.numeric(rep(1,12)))
colnames(dci_sig) <- c("pv")
                     
for (M in 1:12) {
  #subset to sequential months
  k1 = which(ci_daily$month %in% M)
  if (M < 12) {
    k2 = which(ci_daily$month %in% (M+1))
  } else {
    k2 = which(ci_daily$month %in% (M-11))
  }
  # first, compare CI
  s1 <- ci_daily$CI[k1]
  s2 <- ci_daily$CI[k2]
  # remove extreme outliers
  # q1 <- quantile(s1, na.rm = TRUE)
  # iqr1 <- IQR(s1, na.rm = TRUE)
  # s1e <- subset(s1, s1 > (q1[2] - 1.5*iqr1) & s1 < (q1[4] +1.5*iqr1))
  # s2 <- ci_daily$DCI[k2]
  # q2 <- quantile(s2, na.rm = TRUE)
  # iqr2 <- IQR(s2, na.rm = TRUE)
  # s2e <- subset(s2, s2 > (q2[2] - 1.5*iqr2) & s2 < (q2[4] +1.5*iqr2))
  # run t-test
  res <- t.test(s1,s2,alternative = "two.sided", var.equal = FALSE)
  ci_sig$pv[M] <- res$p.value
  
  # now compare delta_CI
  s1 <- ci_daily$DCI[k1]
  s2 <- ci_daily$DCI[k2]
  # remove extreme outliers
  # q1 <- quantile(s1, na.rm = TRUE)
  # iqr1 <- IQR(s1, na.rm = TRUE)
  # s1e <- subset(s1, s1 > (q1[2] - 1.5*iqr1) & s1 < (q1[4] +1.5*iqr1))
  # s2 <- ci_daily$DCI[k2]
  # q2 <- quantile(s2, na.rm = TRUE)
  # iqr2 <- IQR(s2, na.rm = TRUE)
  # s2e <- subset(s2, s2 > (q2[2] - 1.5*iqr2) & s2 < (q2[4] +1.5*iqr2))
  #run t-test
  res <- t.test(s1,s2,alternative = "two.sided", var.equal = FALSE)
  dci_sig$pv[M] <- res$p.value  
}

pvs <- data.frame(cbind(ci_sig$pv,dci_sig$pv))
write_csv(pvs,"ci_pvs.csv")
