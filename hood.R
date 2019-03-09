rm(list= ls())
library(LeeCarteR)

load("E:/Source/LeeCarteR/data/CrudeData.rdata")


p <- data.ct$Male$pm[, 1:85]
x <- data.ct$Male$dea[, 1:85]
fit <- fit_lcm(x, p, link="log")

identify_kt(fit)

rd <- residuals(fit)


obj <- fit

mc <- bootstrap(fit, spe="i1", n_forward=33)





