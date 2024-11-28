sempl_ext <- function(BNaught, OmegaNaught, BTilde, OmegaTilde,
                      y, B, Omega, aelWeight = 0, innerTol = 1e-5, maxInnerIter = 100, maxSearch = 1000, silent = T, gammaFactor = NULL) {

 paramNaught <- matrix(c(BNaught[which(B == 1)],
                                  OmegaNaught[which(Omega== 1 & lower.tri(Omega, diag = T))]),
                                ncol = 1)
 paramTilde <-  init <- matrix(c(BTilde[which(B == 1)],
                                  OmegaTilde[which(Omega== 1 & lower.tri(Omega, diag = T))]),
                                ncol = 1)
  if(is.null(gammaFactor)){
    gammaFactor <- 2 * dim(y)[1]
  }


 theta <- paramTilde + .5 * (paramNaught - paramTilde)
 stepInc <- .25 * (paramNaught - paramTilde)


 theta_lrt <- empLikSem::naive_sempl_lrt(weights = theta,
                         y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges,
                         innerTol = 1e-5, maxInnerIter = 100, aelWeight = aelWeight)
 h_theta <- paramTilde + (1 + theta_lrt / gammaFactor) *  (theta - paramTilde)


counter <- 0
scalingOld <- (h_theta - paramTilde)[1] / (paramNaught - paramTilde)[1]

 while(norm(h_theta - paramNaught) > 1e-2 & counter < maxSearch) {
   counter <- counter + 1

   scalingNew <- (h_theta - paramTilde)[1] / (paramNaught - paramTilde)[1]
   if(sign(log(scalingOld)) != sign(log(scalingNew))){
     stepInc <- .5 * stepInc
   }

   if(scalingNew > 1){
     theta <- theta - stepInc
   } else {
     theta <- theta + stepInc
   }

   theta_lrt <- empLikSem::naive_sempl_lrt(weights = theta,
                                           y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges,
                                           innerTol = 1e-5, maxInnerIter = 100, aelWeight = aelWeight)
   h_theta <- paramTilde + (1 + theta_lrt / gammaFactor) *  (theta - paramTilde)

   scalingOld <- scalingNew

   if(silent == F){
     print(paste("Search: ", counter, "; norm: ", round(norm(h_theta - paramNaught), 2), "; LRT: ",
                 round(theta_lrt,2), "; Scale: ", round(scalingOld, 3), sep = ""))
   }
 }
theta_lrt <- empLikSem::naive_sempl_lrt(weights = theta,
                                        y = out$Y, b_r = out$directEdges, o_r = out$bidirectEdges,
                                        innerTol = 1e-5, maxInnerIter = 100, aelWeight = aelWeight)

return(list(theta = theta, lrt = theta_lrt))
}

