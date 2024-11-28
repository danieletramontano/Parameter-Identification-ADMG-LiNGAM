#' Confidence sets for Causal Discovery
#'
#' Estimate parameters of SEM corresponding to ADMG
#'
#' @param Y a matrix with p columns (corresponding to variables) and each row is an observation. Should be mean zero
#' @param B a p x p binary matrix where B[i,j] = 1 indicates a directed edge j -> i
#' @param O a p x p binary matrix where O[i,j] = O[j,i] = 1 indicates the bidirected edge i <-> j
#' @param innerTol convergence criteria for inner optimization in calculating EL LRT at specific point
#' @param maxInnerIter max iterations for inner optimization in calculating EL LRT
#' @param outerTol convergence criteria for outer optimization searching over parameters
#' @param maxOuterIter max iterations for outer optimization searching over parameters
#' @param naive if FALSE (default) profiles out bidirected parameters. If TRUE, search over all directed and bidirected parameters
#' @param Binit optional initialization for directed edges
#' @param Oinit optional initialization for bidirected edged
#' @param OmegaScale parameter which governs initialization
#' @param aelWEight weight given to extra augmented EL point. Default is 0 (i.e., no AEL point)
#'
#' @return
#' ci is the interval
#' length is the length
sempl <- function(Y, B, O, innerTol = 1e-4, maxInnerIter = 100,
                  outerTol = 1e-5, maxOuterIter = 1000, naive = F, Binit = NULL,
                  Oinit = NULL, OmegaScale = .9, aelWeight = 0){

  ## Get initialization for B if necessary
  if(is.null(Binit)){
    Binit <- matrix(0, dim(Y)[2], dim(Y)[2])
    res <- Y
    for(i in 1:dim(Y)[2]){

      i_parents <- which(B[i, ] == 1)

      if(length(i_parents) > 0){
        mod <- RcppArmadillo::fastLm(Y[, i_parents],  y = Y[, i])
        res[, i] <- mod$residual
        Binit[i, i_parents] <- mod$coefficients
      } else {
        res[, i] <- Y[, i]
      }
    }
  }

  ## Get initialization for Omega if necessary
  if(is.null(Oinit) & naive){
    Oinit <- cov(res) * O
    if(any(eigen(Oinit, only.values = T)$values < 0)){

      ## Check diagonal dominance for each row
      for(i in 1:dim(Y)[2]){
        row_i_sum <- sum(abs(Oinit[i, -i]))

        ## If not diagonally dominant, then scale down
        if(row_i_sum > Om[i, i]){
          Oinit[i, -i] <- Oinit[i, -i] * (Oinit[i, i] / row_i_sum) * OmegaScale
          Oinit[-i, i] <- Oinit[i, -i] # to preserve symmetry
        }
      }
    }
  }


  if(!naive){

    init <- c(Binit[which(B == 1)])
    bfgs_out <- lbfgs::lbfgs(sempl_lrt, sempl_grad, vars = init, y = Y, b_r = B, o_r = O, innerTol = innerTol, maxInnerIter = maxInnerIter, aelWeight = aelWeight,
                             invisible = 1, max_iterations = maxOuterIter, epsilon = outerTol)

    final.mod <- semplC(bfgs_out$par, y = Y, b_r = B, o_r = O, innerTol = innerTol, maxInnerIter = maxInnerIter, aelWeight = 0)

  } else {

    init <- c(Binit[which(B == 1)], Oinit[which(O == 1 & lower.tri(O, diag = T))])
    bfgs_out <- lbfgs::lbfgs(naive_sempl_lrt, naive_sempl_grad, vars = init, y = Y, b_r = B, o_r = O, innerTol = innerTol, maxInnerIter = maxInnerIter, aelWeight = aelWeight, invisible = 1,
                             max_iterations = maxOuterIter, epsilon = outerTol)
    final.mod <- naive_semplC(bfgs_out$par, y = Y, b_r = B, o_r = O, innerTol = innerTol, maxInnerIter = maxInnerIter, aelWeight = 0)
    Binit[which(B == 1)] <- bfgs_out$par[1:sum(B)]
    Oinit[which(O == 1 & lower.tri(O, diag = T))] <- bfgs_out$par[(sum(B) + 1): length(bfgs_out$par)]
  }

  # sigma <- solve(diag(rep(1, dim(Y)[2])) - final.mod$B, final.mod$Omega) %*% t(solve(diag(rep(1, dim(Y)[2])) - final.mod$B))


  return(list(conv = bfgs_out$convergence == 0 & abs(sum(final.mod$weights) - 1) < 1e-5,
              B_hat = final.mod$B,
              O_hat = round(final.mod$Omega, 12),
              weights = final.mod$weights,
              lrt = final.mod$lrt,
              df = final.mod$df,
              dual = final.mod$dual))
              # S_hat = sigma))


}


