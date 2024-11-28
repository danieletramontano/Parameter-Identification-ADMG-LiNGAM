#### Testing Empirical Likelihood for Structural Equation Models ####
library(microbenchmark)
library(lavaan)
library(empLikSem)
library(BCD)
source("run_lavaan.R")




p <- 8
d <- 10
b <- 6


dInd <- 2
args=commandArgs(TRUE)
eval(parse(text=args[[1]]))
print(dInd)

set.seed(1234)
dist.list <- c("gauss", "t", "ln", "gamma")
n.list <- c(100, 250, 500 , 1000)
param.run <- expand.grid(dist.list, n.list)
sim.size <- 500

dist <- as.character(param.run[dInd, 1])
n <- param.run[dInd, 2]

print(dist)
print(n)

sq.err.b <- sq.err.o <- sq.err.sig <-  conv <-  matrix(0, ncol = 6, nrow = sim.size)
p.val <- matrix(0, ncol = 2, nrow = sim.size)

  for(s in 1:sim.size){
    print(s)
    out <- empLikSem::rBAP(p, n, dist = dist, d, b, t.df = 4, cycle = 4)

      out_prof <- NULL
      try(out_prof <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges), silent = T)
      if(is.null(out_prof)){
        out_prof <- list(B.hat = 0, Omega.hat = 0, converged = 4, lr = 1e10)
        sigma_prof <- 0
      } else {
        sigma_prof <- solve(diag(rep(1, p)) - out_prof$B_hat) %*% out_prof$O_hat %*% t(solve(diag(rep(1, p)) - out_prof$B_hat))
      }

      out_adj <- NULL
      try(out_adj <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, aelWeight = log(n)/2), silent = T)
      if(is.null(out_adj)){
        out_adj <- list(B.hat = 0, Omega.hat = 0, converged = 4, lr = 1e10)
        sigma_adj <- 0
      } else {
        sigma_adj <- solve(diag(rep(1, p)) - out_adj$B_hat) %*% out_adj$O_hat %*% t(solve(diag(rep(1, p)) - out_adj$B_hat))

      }

      out_lavaan_gls <- NULL
      try(out_lavaan_gls <- run.lavaan(out$Y, out$directEdges, out$bidirectEdges, type = "GLS"), silent = T)

      if(is.null(out_lavaan_gls)){
        out_lavaan_gls <- list(B.hat = 0, Omega.hat = 0, converged = 0, p.val = 1e-8)
        sigma_gls <- 0
      } else {
        sigma_gls <- solve(diag(rep(1, p)) - out_lavaan_gls$B.hat) %*% out_lavaan_gls$Omega.hat %*% t(solve(diag(rep(1, p)) - out_lavaan_gls$B.hat))
      }

      out_lavaan_wls <- NULL
      try(out_lavaan_wls <- run.lavaan(out$Y, out$directEdges, out$bidirectEdges, type = "WLS"), silent = T)

      if(is.null(out_lavaan_wls)){
        out_lavaan_wls <- list(B.hat = 0, Omega.hat = 0, converged = 0, p.val = 1e-8)
        sigma_wls <- 0
      } else {
        sigma_wls <- solve(diag(rep(1, p)) - out_lavaan_wls$B.hat) %*% out_lavaan_wls$Omega.hat %*% t(solve(diag(rep(1, p)) - out_lavaan_wls$B.hat))
      }

      out_ricf <- ricf(out$directEdges, out$bidirectEdges, t(out$Y), maxIter = 1000, msgs = F)



      hybrid <- semplC(out_ricf$BHat[which(out$directEdges == 1)], out$Y, out$directEdges, out$bidirectEdges, 1e-4, 100, 0)
      if((abs(sum(hybrid$weights) - 1)) < 1e-4){
        sigma_hybrid <- solve(diag(rep(1, p)) - out_ricf$BHat) %*% hybrid$O %*% t(solve(diag(rep(1, p)) - out_ricf$BHat))
        conv_hybrid <- 1
      } else {
        sigma_hybrid <- 0
        conv_hybrid <- 0
      }


      sq.err.b[s, ] <- c(sum(c(out_lavaan_gls$B.hat - out$B)^2),
                         sum(c(out_lavaan_wls$B.hat - out$B)^2),
                         sum(c(out_prof$B_hat - out$B)^2),
                         sum(c(out_adj$B_hat - out$B)^2),
                         sum(c(out_ricf$BHat - out$B)^2),
                         sum(c(out_ricf$BHat - out$B)^2)) / sum(out$B^2)


      sq.err.o[s, ] <- c(sum(c(out_lavaan_gls$Omega.hat[lower.tri(out$Omega, diag = T)] - out$Omega[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(out_lavaan_wls$Omega.hat[lower.tri(out$Omega, diag = T)] - out$Omega[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(out_prof$O_hat[lower.tri(out$Omega, diag = T)] - out$Omega[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(out_adj$O_hat[lower.tri(out$Omega, diag = T)] - out$Omega[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(out_ricf$OmegaHat[lower.tri(out$Omega, diag = T)] - out$Omega[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(hybrid$Omega[lower.tri(out$Omega, diag = T)] - out$Omega[lower.tri(out$Omega, diag = T)])^2)) / sum(out$Omega[lower.tri(out$Omega, diag = T)]^2)






      sq.err.sig[s, ] <- c(sum(c(sigma_gls[lower.tri(out$Omega, diag = T)] - out$Sigma[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(sigma_wls[lower.tri(out$Omega, diag = T)] - out$Sigma[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(sigma_prof[lower.tri(out$Omega, diag = T)] - out$Sigma[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(sigma_adj[lower.tri(out$Omega, diag = T)] - out$Sigma[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(out_ricf$SigmaHat[lower.tri(out$Omega, diag = T)] - out$Sigma[lower.tri(out$Omega, diag = T)])^2),
                         sum(c(sigma_hybrid[lower.tri(out$Omega, diag = T)] - out$Sigma[lower.tri(out$Omega, diag = T)])^2)) / sum(out$Sigma[lower.tri(out$Omega, diag = T)]^2)

      conv[s, ] <- c(out_lavaan_gls$converged & all(eigen(out_lavaan_gls$Omega.hat, only.values = T)$values > 0),
                     out_lavaan_wls$converged & all(eigen(out_lavaan_wls$Omega.hat, only.values = T)$values > 0),
                     (out_prof$conv == 0) & (abs(sum(out_prof$weights) - 1) < 1e-6),
                     (out_adj$conv == 0) & (abs(sum(out_adj$weights) - 1) < 1e-6),
                     out_ricf$Converged,
                     conv_hybrid == 1)
      p.val[s, ] <- c(out_lavaan_gls$p.val, out_lavaan_wls$p.val)
  }
  all.conv <- which(apply(conv, MAR = 1, prod) == 1)
  sq.err.b.total <- colMeans(sq.err.b[all.conv, ])
  sq.err.o.total <- colMeans(sq.err.o[all.conv, ])
  sq.err.sig.total <- colMeans(sq.err.sig[all.conv, ])
  sq.err.b.total.med <-  apply(sq.err.b[all.conv, ], MAR = 2, FUN = median)
  sq.err.o.total.med <- apply(sq.err.o[all.conv, ], MAR = 2, FUN = median)
  sq.err.sig.total.med <- apply(sq.err.sig[all.conv, ], MAR = 2, FUN = median)

  conv.total <- colMeans(conv)

  # write.csv(sq.err.b.total, paste(dist, "_errB_", n, ".csv", sep = ""))
  # write.csv(sq.err.o.total, paste(dist, "_errO_", n, ".csv", sep = ""))
  # write.csv(sq.err.sig.total, paste(dist, "_errS_", n, ".csv", sep = ""))
  #
  # write.csv(sq.err.b.total.med, paste(dist, "_errBmed_", n, ".csv", sep = ""))
  # write.csv(sq.err.o.total.med, paste(dist, "_errOmed_", n, ".csv", sep = ""))
  # write.csv(sq.err.sig.total.med, paste(dist, "_errSmed_", n, ".csv", sep = ""))
  #
  # write.csv(conv.total, paste(dist, "_conv_", n, ".csv", sep = ""))


png("boxplots.png", width= 800, height = 600)
boxplot(sq.err.sig[all.conv, ], names = c("G", "W", "E", "A","R", "H"), main = "T7 with 100 obs")
dev.off()

sq.err.sig.total
