library(empLikSem)
library(lbfgs)

p <- 8
d <- 10
b <- 6


dInd <- 1
args=(commandArgs(TRUE))

for (i in 1:length(args)) {
  eval (parse (text = args[[i]] ))
}
print(dInd)

set.seed(1234)

sim.size <- 500
dist.list <- c("gauss", "t", "ln", "gamma")
n.list <- c(50, 100, 250, 500)

param.run <- expand.grid(dist.list, n.list)
dist <- as.character(param.run[dInd, 1])
n <- param.run[dInd, 2]

print(dist)
print(n)


speed <- sq.err <- conv <- matrix(0, nrow = sim.size, ncol = 6)


    for(s in 1:sim.size){
      print(s)
      out <- rBAP(p, n, dist = dist, d, b, t.df = 7)

      # Profiled Model
      time1 <- microbenchmark::microbenchmark(fit.prof <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, innerTol = 1e-3, outerTol = 1e-4,
                                                        maxInnerIter = 50, maxOuterIter = 500), times = 1)

      # Naive
      time2 <- microbenchmark::microbenchmark(fit.naive <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, naive = T, innerTol = 1e-3, outerTol = 1e-4,
                                                         maxInnerIter = 50, maxOuterIter = 500), times = 1)

      # Adjusted Model
      time3 <- microbenchmark::microbenchmark(fit.prof.adj <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, aelWeight = log(n)/2, innerTol = 1e-3, outerTol = 1e-4,
                                                            maxInnerIter = 50, maxOuterIter = 500), times = 1)
      time4 <- microbenchmark::microbenchmark(fit.naive.adj <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, naive = T, aelWeight = log(n)/2, innerTol = 1e-3, outerTol = 1e-4,
                                                             maxInnerIter = 50, maxOuterIter = 500), times = 1)

      # Hybrid Model
      time5 <- microbenchmark::microbenchmark(fit.prof.hybrid <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, Binit = fit.prof.adj$B_hat, Oinit = fit.prof.adj$O_hat, naive = F, innerTol = 1e-3, outerTol = 1e-4,
                                                               maxInnerIter = 50, maxOuterIter = 500), times = 1)
      time6 <- microbenchmark::microbenchmark(fit.naive.hybrid <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, Binit = fit.naive.adj$B_hat, Oinit = fit.naive.adj$O_hat, naive = T, innerTol = 1e-3, outerTol = 1e-4,
                                                                maxInnerIter = 50, maxOuterIter = 500), times = 1)


      speed[s, ] <- c(time1$time, time2$time, time3$time, time4$time, time5$time, time6$time) / 1e9
      sq.err[s, ] <- c(sum(c(fit.prof$S_hat[lower.tri(out$Sigma)] - out$Sigma[lower.tri(out$Sigma)])^2),
                       sum(c(fit.naive$S_hat[lower.tri(out$Sigma)] - out$Sigma[lower.tri(out$Sigma)])^2),
                       sum(c(fit.prof.adj$S_hat[lower.tri(out$Sigma)] - out$Sigma[lower.tri(out$Sigma)])^2),
                       sum(c(fit.naive.adj$S_hat[lower.tri(out$Sigma)] - out$Sigma[lower.tri(out$Sigma)])^2),
                       sum(c(fit.prof.hybrid$S_hat[lower.tri(out$Sigma)] - out$Sigma[lower.tri(out$Sigma)])^2),
                       sum(c(fit.naive.hybrid$S_hat[lower.tri(out$Sigma)] - out$Sigma[lower.tri(out$Sigma)])^2)) / sum(out$Sigma[lower.tri(out$Sigma)]^2)
      conv[s, ] <- c(fit.prof$conv, fit.naive$conv, fit.prof.adj$conv, fit.naive.adj$conv, fit.prof.hybrid$conv, fit.naive.hybrid$conv)
    }
    all.conv <- which(apply(conv, MAR = 1, prod) == 1)

    if(length(all.conv) > 0){
      speed.total <- colMeans(speed[all.conv, , drop = F])
      sq.err.total <- colMeans(sq.err[all.conv, , drop = F])
    } else {
      speed.total <- rep(0, 6)
      sq.err.total <- rep(0, 6)
    }
    conv.total <- colMeans(conv)

    # write.csv(speed.total, paste(dist, "_speed", n, ".csv", sep = ""))
    # write.csv(sq.err.total, paste(dist, "_sqErr", n, ".csv", sep = ""))
    # write.csv(conv.total, paste(dist, "_conv", n, ".csv", sep = ""))
    # write.csv(agree.total, paste(dist, "_agree", n, ".csv", sep = ""))




