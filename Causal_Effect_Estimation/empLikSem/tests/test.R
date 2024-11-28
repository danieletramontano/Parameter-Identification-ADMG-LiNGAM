library(empLikSem)
library(lbfgs)

p <- 8
n <- 50
d <- 12
b <- 6

rec <- matrix(0, nrow = sim.size, ncol = 2)
out <- rBAP(p, n, dist = "gauss", d, b)



init <- out$B[which(out$directEdges == 1)]
outGrad <- numDeriv::grad(sempl_lrt, init, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, innerTol = 1e-4,
                          maxInnerIter =  500, aelWeight = log(n)/2)
outGrad1 <- sempl_grad(init, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, innerTol = 1e-4, maxInnerIter =  100, aelWeight = log(n)/2)
outGrad / outGrad1





bfgs_out <- lbfgs(sempl_lrt, sempl_grad, vars = init, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, innerTol = 1e-4, maxInnerIter =  100, aelWeight = log(n)/2)

mod.fit.adj <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, aelWeight = 0)
mod.fit.ael <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, aelWeight = log(n)/2)

mod.fit.ael$lrt - mod.fit$lrt
mod.fit$B_hat - mod.fit.ael$B_hat


varEL <- empLikSem::var.el(Y = out$Y, B = out$directEdges, Omega = out$bidirectEdges, B.hat = mod.fit$B_hat, Omega.hat = mod.fit$O_hat, weights = c(mod.fit$weights))

(est - c(init, out$O[out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T)]))

/ diag(varEL)




weights <- c(out$B[which(out$directEdges == 1)], out$Omega[which(out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T))])
dg <- naive_sempl_dg(weights, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, innerTol = 1e-4, maxInnerIter =  100)
mod.fit <- naive_semplC(weights, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, 1e-4, 100)
sum(log(mod.fit$weights))





num_deriv <- numDeriv::grad(naive_sempl_lrt, weights, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, innerTol = 1e-4, maxInnerIter =  100)
round(num_deriv - naive_grad, 5)



bfgs_out.naive <- lbfgs(naive_sempl_lrt, naive_sempl_grad, vars = weights, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, innerTol = 1e-4, maxInnerIter =  100)
final.mod.naive <- naive_semplC(bfgs_out.naive$par, y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges, innerTol = 1e-4, maxInnerIter =  100)




mod.fit.naive <- sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges, naive = T)



round(mod.fit.naive$O_hat - mod.fit$O_hat, 5)
round(mod.fit.naive$B_hat - mod.fit$B_hat, 5)

sim.size <- 1000
rec <- matrix(0, nrow = sim.size, ncol = 4)

p <- 8
n <- 500
d <- 10
b <- 6
for(i in 1:sim.size){
  cat(i)
  out <- rBAP(p, n, dist = "gamma", d, b)
  init <- matrix(c(out$B[which(out$directEdges == 1)],
                   out$Omega[which(out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T))]),
                 ncol = 1)
  mele <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges)
  adj.out <- empLikSem::naive_semplC(weights = init,
                                     y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges,
                                     innerTol = 1e-5, maxInnerIter = 100, aelWeight = .5 * log(n))
  unadj.out <- empLikSem::naive_semplC(weights = init,
                                       y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges,
                                       innerTol = 1e-5, maxInnerIter = 100, aelWeight = 0)
  extend.out <- sempl_ext(out$B, out$Omega, mele$B_hat, mele$O_hat,
                          out$Y, out$directEdges, out$bidirectEdges,aelWeight = 0.0, maxSearch = 1000, silent =T)
  extend.out2 <- sempl_ext(out$B, out$Omega, mele$B_hat, mele$O_hat,
                          out$Y, out$directEdges, out$bidirectEdges,aelWeight = 0.0, maxSearch = 1000, gammaFactor = n, silent =T)
  extend.out.adj <- sempl_ext(out$B, out$Omega, mele$B_hat, mele$O_hat,
                              out$Y, out$directEdges, out$bidirectEdges,aelWeight = log(n)/2, maxSearch = 1000, silent =T)



  rec[i, ] <- c(extend.out$lrt, extend.out.adj$lrt, adj.out$lrt, unadj.out$lrt) - mele$lrt
}


colMeans(rec < qchisq(.95, df = b + d+ p))
