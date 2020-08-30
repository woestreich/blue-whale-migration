library(tidyverse)

rm(list = ls())

#2018
lunges <- read.csv("lungesTDR2018.csv",header = FALSE)
colnames(lunges) <- c("day","night")
lunges <- lunges[c(1:10),]

calls <- read.csv("callsTDR2018.csv",header = FALSE)
colnames(calls) <- c("day","night")
calls_for <- calls[c(1:10),]
calls_migr <- calls[c(11:31),]

s1 <- calls_for$night - calls_for$day
s2 <- calls_migr$night - calls_migr$day

res <- t.test(s1,s2,alternative = "two.sided", var.equal = FALSE)

#2019
rm(list = ls())

lunges <- read.csv("lungesTDR2019.csv",header = FALSE)
colnames(lunges) <- c("day","night")
lunges <- lunges[c(1:7),]

calls <- read.csv("callsTDR2019.csv",header = FALSE)
colnames(calls) <- c("day","night")
calls_for <- calls[c(1:7),]
calls_migr <- calls[c(8:17),]

s1 <- calls_for$night - calls_for$day
s2 <- calls_migr$night - calls_migr$day

res <- t.test(s1,s2,alternative = "two.sided", var.equal = FALSE)
