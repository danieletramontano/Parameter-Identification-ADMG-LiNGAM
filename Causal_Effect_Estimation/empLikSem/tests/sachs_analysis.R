## Add bidirected Edge

rm(list = ls())
data <- read.csv("data/sachs_data_cd3cd28.csv")


Y <- matrix(t(data), ncol = 11)

colnames(Y) <- names(data)

subY <- log(Y)
subYscale <- scale(subY)



B <- matrix(0, nrow = dim(subY)[2], ncol = dim(subY)[2])
rownames(B) <- colnames(B) <- colnames(subY)
B["PKC", "plcg"] <- 1
B["PKC", "PIP2"] <- 1
B["PIP2", "plcg"] <- 1
B["plcg", "PIP3"] <- 1
B["pakts473", "PIP3"] <- 1
B["pjnk", "PKC"] <- 1
B["praf", "PKC"] <- 1
B["P38", "PKC"] <- 1
B["pakts473", "PKA"] <- 1
B["pjnk", "PKA"] <- 1
B["P38", "PKA"] <- 1
B["praf", "PKA"] <- 1
B["p44.42", "PKA"] <- 1
B["p44.42", "pmek"] <- 1
B["pmek", "praf"] <- 1




O <- matrix(0, nrow = dim(subY)[2], ncol = dim(subY)[2])
rownames(O) <- colnames(O) <- colnames(subY)
O["PIP2", "PIP3"] <- 1
O["PIP2", "praf"] <- 1
O["PIP3", "praf"] <- 1
O <- O + t(O) + diag(rep(1, dim(subY)[2]))

O.submod <- O
O["pakts473", "PIP2"] <- O["PIP2", "pakts473"] <-1

B.submod <- B


out.ricf <- BCD::ricf(B, O,t(subYscale))
out.sempl <- empLikSem::sempl(subYscale, B, O)

out.ricf.submod <- BCD::ricf(B.submod, O.submod, t(subYscale))
out.sempl.submod <- empLikSem::sempl(subYscale, B.submod, O.submod)


sempl.pval <- pchisq((out.sempl.submod$lrt - out.sempl$lrt), df = 1,lower.tail = F)

ricf.pval <- pchisq(2* (sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf$S, log = T)) - sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf.submod$S, log = T))),
                    df = 1, lower.tail = F)

(out.sempl.submod$lrt - out.sempl$lrt)
sempl.pval
2* (sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf$S, log = T)) - sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf.submod$S, log = T)))
ricf.pval



# colnames(B)
#
# dual.ranked <- which(lower.tri(O.submod) & O.submod != 1, arr.ind = T)[order(abs(out.sempl.submod$dual[-c(1:dim(B)[1])]), decreasing = F), ]
cbind(colnames(B)[dual.ranked[, 1]], colnames(B)[dual.ranked[, 2]])
################################################
## Directed Cycle

rm(list = ls())

data <- read.csv("data/sachs_data_cd3cd28.csv")



Y <- matrix(t(data), ncol = 11)

colnames(Y) <- names(data)

subY <- log(Y)
subYscale <- scale(subY)


B <- matrix(0, nrow = dim(subY)[2], ncol = dim(subY)[2])
rownames(B) <- colnames(B) <- colnames(subY)
B["PKC", "plcg"] <- 1
B["PKC", "PIP2"] <- 1
B["PIP2", "plcg"] <- 1
B["plcg", "PIP3"] <- 1
B["pakts473", "PIP3"] <- 1
B["pjnk", "PKC"] <- 1
B["praf", "PKC"] <- 1
B["P38", "PKC"] <- 1
B["pakts473", "PKA"] <- 1
B["pjnk", "PKA"] <- 1
B["P38", "PKA"] <- 1
B["praf", "PKA"] <- 1
B["p44.42", "PKA"] <- 1
B["p44.42", "pmek"] <- 1
B["pmek", "praf"] <- 1





O <- matrix(0, nrow = dim(subY)[2], ncol = dim(subY)[2])
rownames(O) <- colnames(O) <- colnames(subY)
O["PIP2", "PIP3"] <- 1
O["PIP2", "praf"] <- 1
O["PIP3", "praf"] <- 1
O <- O + t(O) + diag(rep(1, dim(subY)[2]))

O.submod <- O

B.submod <- B
B["PKA", "pmek"] <- 1

colnames(B)

out.ricf <- BCD::ricf(B, O,t(subYscale))
out.sempl <- empLikSem::sempl(subYscale, B, O)

out.ricf.submod <- BCD::ricf(B.submod, O.submod, t(subYscale))
out.sempl.submod <- empLikSem::sempl(subYscale, B.submod, O.submod)


sempl.pval <- pchisq((out.sempl.submod$lrt - out.sempl$lrt), df = 1,lower.tail = F)
ricf.pval <- pchisq(2* (sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf$S, log = T)) - sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf.submod$S, log = T))),
                    df = 1, lower.tail = F)


(out.sempl.submod$lrt - out.sempl$lrt)
2* (sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf$S, log = T)) - sum(mvtnorm::dmvnorm(subYscale, sigma = out.ricf.submod$S, log = T)))

sempl.pval
ricf.pval

