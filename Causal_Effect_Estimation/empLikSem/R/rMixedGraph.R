rBAP <- function(p, n, dist, d, b, t.df = 4, cycle = 0){

  ### Build check for d + b too large ###

  B <- Omega <- matrix(0, nrow = p, ncol = p)


  if(cycle > 0){
    B[matrix(c(2:cycle, 1:(cycle - 1)), ncol = 2)] <- 1
    B[1, cycle] <- 1
  }

  # Draw edges from uniform(.2, 1)
  B[sample(which(lower.tri(B) & B != 1), size = d)] <- 1
  B.true <-  sample(c(-1, 1), size = p^2, replace = T) * runif(p^2, .2, 1) * B

  #
  Omega[sample(which(lower.tri(Omega) & B != 1), size = b)] <- 1

  Omega.true <- sample(c(-1, 1), size = p^2, replace = T) * runif(p^2, .3, .8) * Omega
  Omega.true <- Omega.true + t(Omega.true)
  diag(Omega.true) <- 1 + rowSums(abs(Omega.true)) + rexp(p)

  Omega <- Omega + t(Omega) + diag(rep(1, p))


  if(dist == "gauss"){
    errs <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Omega.true)
  } else if(dist == "t") {
    errs <- mvtnorm::rmvt(n, sigma = Omega.true * (t.df - 2)/ (t.df), df = t.df)
  } else if(dist == "ln_t") {
    omega.cor <- diag(1/sqrt(diag(Omega.true))) %*% Omega.true %*% diag(1/sqrt(diag(Omega.true)))
    Omega.true <- exp(1) * (exp(omega.cor) - 1)
    errs <- mvtnorm::rmvt(n, sigma = Omega.true * (t.df - 2)/ (t.df), df = t.df)
  } else if(dist == "ln") {
    omega.cor <- diag(1/sqrt(diag(Omega.true))) %*% Omega.true %*% diag(1/sqrt(diag(Omega.true)))
    Z <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = omega.cor)
    Omega.true <- exp(1) * (exp(omega.cor) - 1)
    errs <- (exp(Z) - exp(1/2))
  } else if (dist == "gamma"){

    ### Variance Terms ###
    gamma.var <- diag(Omega.true) - rowSums(abs(Omega.true - diag(diag(Omega.true))))
    gamma.mean <- gamma.var
    errs <- sapply(gamma.var, FUN = function(x){rgamma(n, shape = x, scale = 1)})

    ### Covariance Terms ###
    non.zeros <- which(Omega == 1 & lower.tri(Omega), arr.ind = T)
    for(i in 1:dim(non.zeros)[1]){
      cov.term <- Omega.true[non.zeros[i, , drop = F]]
      latent.err <- rgamma(n, shape = abs(cov.term), scale = 1)

      # Either both are same sign, or opposite signs based on sign of covariance
      if(cov.term > 0){
        latent.add <- matrix(c(latent.err, latent.err), nrow = n)
        mean.add <- c(abs(cov.term), abs(cov.term))
      } else {
        latent.add <- matrix(c(latent.err, -latent.err), nrow = n)
        mean.add <- c(abs(cov.term), -abs(cov.term))
      }
      # flip direction of skew at random
      flip <- sample(c(1,-1), size = 1)
      latent.add <- flip * latent.add
      mean.add <- flip * mean.add

      gamma.mean[c(non.zeros[i, ])] <-  gamma.mean[c(non.zeros[i, ])] + mean.add
      errs[, non.zeros[i, 1]] <- errs[, non.zeros[i, 1]] + latent.add[, 1]
      errs[, non.zeros[i, 2]] <- errs[, non.zeros[i, 2]] + latent.add[, 2]
    }

    errs <-t(t(errs) - gamma.mean)

  }

  Y <- solve(diag(rep(1, p)) - B.true, t(errs))
  Y <- t(Y)
 sigma <- solve(diag(rep(1, p)) - B.true, Omega.true) %*% t(solve(diag(rep(1,p)) - B.true))

 return(list(B = B.true, Omega = Omega.true, directEdges = B, bidirectEdges = Omega,
        errs = errs, Y = Y, Sigma = sigma))
}

