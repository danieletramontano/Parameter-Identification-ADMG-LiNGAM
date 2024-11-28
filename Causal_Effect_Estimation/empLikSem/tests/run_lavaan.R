library(lavaan)

run.lavaan <- function(Y, B, Omega, type = "ML"){
  mod <- ""
  V <- dim(B)[1]
  for(i in 1:V){

    if(sum(B[i,]) > 0){
      newLinesB <- paste(paste("y",i, sep =""), "~", paste(paste("y", which(B[i,]==1), sep = ""), collapse = "+") )
      mod <- paste(mod, newLinesB, "; ", sep = "")
    }

    newLinesO <- paste(paste("y",i, sep =""), "~~", paste(paste("y", which(Omega[i,1:i]==1), sep = ""), collapse = "+"))
    mod <- paste(mod, newLinesO, "; ", sep = "")

  }

  dat <- data.frame(Y)
  names(dat) <- paste("y", c(1:dim(Y)[2]), sep = "")

  st <- proc.time()
  fit <- lavaan::lavaan(mod, data = dat, fixed.x = F, estimator = type)
  end <- proc.time()
  parOutput <- fit@ParTable


  ## Get row and col of parameters
  rowOut <- sapply(parOutput$lhs, function(x){as.numeric(substring(x, 2))})
  colOut <- sapply(parOutput$rhs, function(x){as.numeric(substring(x, 2))})

  ## which rows of ParOutput$est are for the directed edges
  b.indices <- which(parOutput$op == "~")

  B.hat <- Omega.hat <-matrix(0, nrow = V, ncol = V)
  B.hat[cbind(rowOut[b.indices], colOut[b.indices])] <- parOutput$est[b.indices]

  Omega.hat[cbind(rowOut[-b.indices], colOut[-b.indices])] <- parOutput$est[-b.indices]
  Omega.hat <- Omega.hat + t(Omega.hat) - diag(diag(Omega.hat))

  ## Variance of estimates
  thetaVar <- vcov(fit)


  ### correctly swap row or col of omega indices so that they correspond to a lower triangular matrix
  o.indices <- cbind(rowOut[-b.indices], colOut[-b.indices])
  o.indices <- t(apply(o.indices, MAR = 1,
                     function(x){
                       if(x[1] > x[2]){
                       return(x)
                       } else {
                           return(x[c(2,1)])
                         }
                       }))

  ### Order variance matrix so that it is vec(B) and vec(Omega)
  varOrder <- c(b.indices[order(colOut[b.indices], rowOut[b.indices])],
                c(1:length(parOutput$est))[-b.indices][order(o.indices[, 2], o.indices[, 1])])


  return(list(B.hat = B.hat, Omega.hat = Omega.hat, time = (end-st)[3],
              converged = fit@Fit@converged, p.val = fit@Fit@test[[1]]$pvalue,
              fit = fit, thetaVar = thetaVar[varOrder, varOrder]))
}
#
# out_lavaan <- run.lavaan(Y, B, Omega)
# # out_sempl <- sempl(Y, B, Omega)
#
#
# V <- 5
# n <- 100
# B <- matrix(c(0, 0, 0, 0, 1,
#               1, 0, 0, 0, 0,
#               0, 1, 0, 0, 0,
#               0, 0, 1, 0, 0,
#               0, 0, 0, 1, 0), byrow = T, nrow = 5)
# Omega <- diag(rep(1,V))
#
#
# B.true <- (rbinom(size = 1, n = V^2, prob = .5) *2 - 1) * runif(V^2, .2, 1) * B
# temp <- solve(diag(rep(1,V)) - B.true)
#
# Omega.true <- matrix(0, nrow = V, ncol = V)
# Omega.true[(Omega == 1) & lower.tri(Omega)] <-
# Omega.true <- Omega.true + t(Omega.true)
# diag(Omega.true) <- 5
#
# errs <- matrix(rpois(V*n, 5), nrow = V)
# Y <- temp %*% errs
#
# out_lavaan <- run.lavaan(Y, B, Omega)
# out_sempl <- sempl(Y, B, Omega)
#
# sum(c(out_lavaan$B.hat - B.true)^2)
# sum(c(out_sempl$B.hat - B.true)^2)
# sum(c(out_lavaan$Omega.hat - Omega.true)^2)
# sum(c(out_sempl$Omega.hat - Omega.true)^2)
#
