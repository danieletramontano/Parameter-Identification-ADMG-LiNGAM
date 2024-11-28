#### Testing Empirical Likelihood for Structural Equation Models ####
library(lavaan)
source("run_lavaan.R")
library(empLikSem)
library(lbfgs)
library(BCD)


p <- 6
d <- 8
b <- 4


dInd <- 1
args=(commandArgs(TRUE))

for (i in 1:length(args)) {
  eval (parse (text = args[[i]] ))
}
print(dInd)


dist.list <- c("gauss", "t", "ln", "gamma")
n.list <- c(250, 500, 1000, 2500, 5000)
param.run <- expand.grid(dist.list, n.list)
sim.size <- 500
param.run <- expand.grid(dist.list, n.list)
dist <- as.character(param.run[dInd, 1])
n <- param.run[dInd, 2]


joint.total <- matrix(0, nrow = length(n.list), ncol = 8)
length.o.total <- length.b.total <- matrix(0, nrow = length(n.list), ncol = 5)

set.seed(1111)


  joint <- matrix(0, ncol = 8, nrow = sim.size)
  length.o <- length.b <-  matrix(0, ncol = 5, nrow = sim.size)
  conv <- matrix(0, ncol = 6, nrow = sim.size)
  for(s in 1:sim.size){
    print(s)
    out <- rBAP(p, n, dist = dist, d, b)

    trueParam <- matrix(c(out$B[which(out$directEdges == 1)],
                    out$Omega[which(out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T))]),
                    nrow = b + d + p)


    ### EL Estimates ###
    out_prof <- NULL
    try(out_prof <- empLikSem::sempl(Y = out$Y, B = out$directEdges, O = out$bidirectEdges), silent = T)
    if(is.null(out_prof)){
      out_prof <- list(B.hat = 0, Omega.hat = 0, converged = 4, lr = 1e10)
      elJoint <- 0
      elVar <- matrix(0, nrow = b + p + d, ncol = b + p + d)
    } else {
      elEst <- matrix(c(out_prof$B_hat[which(out$directEdges == 1)],
                 out_prof$O_hat[which(out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T))]),
                 nrow = b + d+ p)
      elVar <- empLikSem::var.el(Y = out$Y, B = out$directEdges, Omega = out$bidirectEdges,
                                 B.hat = out_prof$B_hat, Omega.hat = out_prof$O_hat, weights = c(out_prof$weights))
      elJoint <- t(elEst - trueParam) %*% solve(elVar$thetaVar) %*% (elEst - trueParam) < qchisq(.95, df = p + d + b)
    }


    ### EL LRT ###
    out_el_direct <- NULL
    try(out_el_direct <- empLikSem::naive_semplC(weights = trueParam,
                                           y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges,
                                           innerTol = 1e-5, maxInnerIter = 100, aelWeight = 0), silent = T)
    if(is.null(out_el_direct)){
      elDirect <- 0
    } else {
      elDirect <- out_el_direct$lrt - out_prof$lrt < qchisq(.95, df = b + d + p)
    }


    ### EL LRT ###
    out_el_adj <- NULL
    try(out_el_adj <- empLikSem::naive_semplC(weights = trueParam,
                                                 y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges,
                                                 innerTol = 1e-5, maxInnerIter = 100, aelWeight = log(n)/2), silent = T)
    if(is.null(out_el_adj)| is.null(out_prof)){
      elAdj <- 0
    } else {
      elAdj <- (out_el_adj$lrt - out_prof$lrt) < qchisq(.95, df = b + d+ p)
    }

    ### EL LRT ###
    out_el_ext <- NULL
    try(out_el_ext <- empLikSem::sempl_ext(BNaught = out$B, OmegaNaught = out$Omega, BTilde = out_prof$B_hat, OmegaTilde = out_prof$O_hat,
                                           y = out$Y, B = out$directEdges, Omega = out$bidirectEdges, aelWeight = 0, innerTol = 1e-4, maxInnerIter = 100, maxSearch = 1000, gammaFactor = 2 * n, silent = T), silent = T)
    if(is.null(out_el_ext) | is.null(out_prof)){
      elExt <- 0
    } else {
      elExt <- (out_el_ext$lrt - out_prof$lrt) < qchisq(.95, df = b + d+ p)
    }


    ### Lavaan ###
    out_lavaan_gls <- NULL
    try(out_lavaan_gls <- run.lavaan(out$Y, out$directEdges, out$bidirectEdges, type = "GLS"), silent = T)

    if(is.null(out_lavaan_gls)){
      out_lavaan_gls <- list(B.hat = 0, Omega.hat = 0, converged = 0, p.val = 1e-8, thetaVar = matrix(0, nrow = b + p + d, ncol = b + p + d))
      glsJoint <- 0
    } else {
      glsEst <- matrix(c(out_lavaan_gls$B.hat[which(out$directEdges == 1)],
                         out_lavaan_gls$Omega.hat[which(out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T))]),
                          nrow = b + d + p)

    glsJoint <- t(glsEst - trueParam) %*% solve(out_lavaan_gls$thetaVar) %*% (glsEst - trueParam) < qchisq(.95, df = p + d+ b)
  }

    out_lavaan_wls <- NULL
    try(out_lavaan_wls <- run.lavaan(out$Y, out$directEdges, out$bidirectEdges, type = "WLS"), silent = T)

    if(is.null(out_lavaan_wls)){
      out_lavaan_wls <- list(B.hat = 0, Omega.hat = 0, converged = 0, p.val = 1e-8, thetaVar = matrix(0, nrow = b + p + d, ncol = b + p + d))
      wlsJoint <- 0
    }  else {
      wlsEst <- matrix(c(out_lavaan_wls$B.hat[which(out$directEdges == 1)],
                         out_lavaan_wls$Omega.hat[which(out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T))]),
                       nrow = b + d+ p)
      wlsJoint <- t(wlsEst - trueParam) %*% solve(out_lavaan_wls$thetaVar) %*% (wlsEst - trueParam) < qchisq(.95, df = p + d+ b)
    }


    ### RICF ###
    out_ricf <- ricf(out$directEdges, out$bidirectEdges, t(out$Y), maxIter = 1000, msgs = F)


    ricfEst <- matrix(c(out_ricf$BHat[which(out$directEdges == 1)],
                           out_ricf$OmegaHat[which(out$bidirectEdges == 1 & lower.tri(out$bidirectEdges, diag = T))]),
                         nrow = b + d+ p)
    expVar <- var.ricf(t(out$Y), B = out$directEdges, Omega = out$bidirectEdges,
             B.hat = out_ricf$BHat, Omega.hat = out_ricf$OmegaHat, type = "expected")
    expJoint <- t(ricfEst - trueParam) %*% solve(expVar) %*% (ricfEst - trueParam) < qchisq(.95, df = p + d+ b)

    sandVar <- var.ricf(t(out$Y), B = out$directEdges, Omega = out$bidirectEdges,
                       B.hat = out_ricf$BHat, Omega.hat = out_ricf$OmegaHat, type = "sandwich")
    sandJoint <- t(ricfEst - trueParam) %*% solve(sandVar) %*% (ricfEst - trueParam) < qchisq(.95, df = p + d+ b)




    joint[s, ] <- c(elJoint, elDirect, elAdj, elExt, glsJoint, wlsJoint, expJoint, sandJoint)

    ci.length <- apply(cbind(diag(elVar$thetaVar), diag(out_lavaan_gls$thetaVar), diag(out_lavaan_wls$thetaVar), diag(expVar), diag(sandVar)), MAR = 1,
                       function(x){
                         sqrt(x) / min(sqrt(x))})
    length.o[s, ] <- rowMeans(ci.length[, (d+1):(b + d + p)])
    length.b[s, ] <- rowMeans(ci.length[, 1:d])

    conv[s, ] <- c(out_lavaan_gls$converged & all(eigen(out_lavaan_gls$Omega.hat, only.values = T)$values > 0),
                   out_lavaan_wls$converged & all(eigen(out_lavaan_wls$Omega.hat, only.values = T)$values > 0),
                   out_prof$conv,
                   out_el_direct$conv,
                   out_el_adj$conv,
                   out_ricf$Converged)
  }

  all.conv <- which(apply(conv, MAR = 1, prod) == 1)
  joint.total <- colMeans(joint[all.conv, ])
  length.o.total <- colMeans(length.o[all.conv, ])
  length.b.total <- colMeans(length.b[all.conv, ])

  # write.csv(joint.total, paste(dist, "_joint", n, ".csv", sep = ""))
  # write.csv(length.o.total, paste(dist, "_lengthO", n, ".csv", sep = ""))
  # write.csv(length.b.total, paste(dist, "_lengthB", n, ".csv", sep = ""))




