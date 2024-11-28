#' Empirical Likelihood Methods for Linear Structural Equation Models
#'
#' Calculate assymptotic variance for non-structural zeros in B and Omega matrix. The ordering of the estimates
#' in the variance matrix follows standard vec(-) convention indexing over rows first then columns. The elements of B
#' are ordered before the elements of Omega. Duplicate elements of Omega are not included.
#'
var.el <- function(Y, B, Omega, B.hat, Omega.hat, weights)
{
  V <- dim(Y)[2]
  n <- dim(Y)[1]
  b_ind <- which(B == 1, arr.ind = T)
  o_ind <- which(Omega == 1 & lower.tri(Omega, diag = T), arr.ind = T)

  errors <- Y - Y %*% t(B.hat)

    # pre-allocate memory (note that dg.i is cleared out in each iteration of the loop)
  dg <- matrix(0, nrow = V + V * (V + 1) /2,
                 ncol =  dim(b_ind)[1] + dim(o_ind)[1])
  g.var <- matrix(0, nrow = V + V * (V + 1) /2,
                    ncol = V + V * (V + 1) /2)


    ## Mean constraints
  for(i in 1:V){
      pa.i <- which(b_ind[, 1] == i)
      if(length(pa.i) > 0){
        dg[i, pa.i] <- -apply(Y[, b_ind[pa.i, , drop = F][, 2], drop = F] * weights, MAR = 2, sum)
      }
    }

    count <- V + 1
    ## Covariance Constraints ##
    for(j in 1:V){
      for(i in j:V){
        if(i != j){

          pa.i <- which(b_ind[, 1] == i)
          if(length(pa.i) > 0){
            dg[count, pa.i] <- -apply(Y[, b_ind[pa.i, ,drop = F][, 2], drop = F] * errors[, j] * weights, MAR = 2, sum)
          }

          pa.j <- which(b_ind[, 1] == j)
          if(length(pa.j) > 0){
            dg[count, pa.j] <- -apply(Y[, b_ind[pa.j, , drop = F][, 2], drop = F] * errors[, i] * weights, MAR = 2, sum)
          }
        } else {
          pa.i <- which(b_ind[, 1] == i)
          if(length(pa.i) > 0){
            dg[count, pa.i] <- -2 * apply(Y[, b_ind[pa.i, , drop = F][, 2], drop = F] * errors[, i] * weights , MAR = 2, sum)
          }
        }

        if(Omega[i, j] == 1){
          dg[count, dim(b_ind)[1] + which(o_ind[, 1] == i & o_ind[, 2] == j)] <- -1
        }

         count <- count + 1

        }

      }

  g <-  cbind(errors, apply(which(lower.tri(Omega, diag = T), arr.ind = T), MAR = 1,
                                function(x){errors[, x[1]] * errors[, x[2]] - Omega.hat[x[1], x[2]]}))
  g.var <- t(g) %*% (g * weights)

  return(list(dg = dg, g = g, gVar = g.var, thetaVar = solve(t(dg) %*% solve(g.var) %*% dg) / n))

}




#' Maximum Likelihood Method for Linear Structural Equation Models
#'
#' Calculate Fisher information for non-structural zeros in B and Omega matrix. The ordering of the estimates
#' in the Information matrix follows standard vec(-) convention indexing over rows first then columns. The elements of B
#' are ordered before the elements of Omega. Duplicate elements of Omega are not included.
#'
#' @param S V by V matrix corresponding to the sample covariance
#' @param B V by V matrix with {0,1} giving structure of directed edges
#' @param Omega V by V matrix with {0,1} giving structure of bi-directed edges
#' @param B.hat V by V matrix giving estimated edges weights for directed edges
#' @param Omega.hat V by V matrix giving estimated edge weights for bi-directed edges
#' @param type string describing which Fisher Information to calculate. Options are "expected" or "observed".
#' @return The inverse (scaled by n) Fisher information matrix as derived by Fox and Drton 2014 or the Huber-White misspecified model covariance estimate
#' @export
var.ricf <- function(Y, B, Omega,  B.hat, Omega.hat, type = "expected")
{

  # Model based SE's with either expected or observed information
  if(type == "expected" | type == "observed"){
    V <- dim(B)[1]
    I <- matrix(0, nrow = sum(B) + (sum(Omega)-V)/2 + V,
                ncol = sum(B) + (sum(Omega)-V)/2 + V)
    n <- dim(Y)[2]
    S <- Y %*% t(Y) / n


    #### Setup permutation matrix P for B
    P <- matrix(0, nrow = V^2, ncol = sum(B))
    temp <- matrix(c(which(c(B)==1),1:sum(B)), nrow = sum(B), ncol = 2)
    P[temp] <- 1


    ### Setup permutation matrix Q for Omega
    Omega.temp <- matrix(0, nrow = V, ncol = V)
    Omega.temp[lower.tri(Omega.temp, diag = T) & Omega == 1] <- c(1:((sum(Omega) - V)/2 + V))
    Omega.temp <- Omega.temp + t(Omega.temp) - diag(diag(Omega.temp))
    Q <- matrix(0, nrow = V^2, ncol = (sum(Omega)-V)/2 + V)
    Q[cbind(c(1:V^2), c(Omega.temp))] <- 1


    ### Setup Transpositon matrix such that vec(t(I-B)) = Tr %*% vec(I-B)
    Tr <- matrix(0, nrow = V^2, ncol = V ^2)
    Tr[matrix(c(1:V^2, (rep(c(1:V), V)-1) *V + rep(c(1:V), each = V)), ncol = 2)] <- 1

    ### precompute quantities of interest
    omega.inv <- solve(Omega.hat)
    eye.minus.b.inv <- solve(diag(rep(1, V)) - B.hat)


    if(type == "expected"){
      ## Expected Information where we assume E(S) = Sigma

      # information for beta elements
      I[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P

      # information for omega elements
      I[-c(1:sum(B)),-c(1:sum(B))] <- 1/2 * t(Q) %*% (omega.inv %x% omega.inv) %*% Q

      # co-information
      I[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ( eye.minus.b.inv %x% omega.inv) %*% Q
      I[-c(1:sum(B)),c(1:sum(B))] <- t(I[c(1:sum(B)),-c(1:sum(B))])

    } else if (type == "observed"){
      ## Observed Information where we do not assume E(S) = Sigma

      # precompute term
      back_term <- omega.inv %*% (diag(rep(1, V)) - B.hat) %*% S %*% t((diag(rep(1, V)) - B.hat)) %*% omega.inv

      # information for beta elements
      I[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P

      # information for omega elements (notice difference between this and expected info above)
      I[-c(1:sum(B)),-c(1:sum(B))] <- -1/2 * t(Q) %*%((omega.inv %x% omega.inv) - (omega.inv %x% back_term) -
                                                        (back_term %x% omega.inv)) %*% Q

      # co-information for omega elements
      I[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ((S %*% t(diag(rep(1, V)) - B.hat) %*% omega.inv) %x% omega.inv) %*% Q
      I[-c(1:sum(B)),c(1:sum(B))] <- t(I[c(1:sum(B)),-c(1:sum(B))])
    }

    # return inverse of FI
    return(solve(I)/n)

  } else if (type == "sandwich") {
    #### Sandwich Variance based on misspecified model ####
    n <- dim(Y)[2]
    V <- dim(B)[1]

    S <- Y %*% t(Y) / dim(Y)[2]

    ## Setup permutation matrix P for B
    P <- matrix(0, nrow = V^2, ncol = sum(B))
    temp <- matrix(c(which(c(B)==1),1:sum(B)), nrow = sum(B), ncol = 2)
    P[temp] <- 1


    ## Setup permutation matrix Q for Omega
    Omega.temp <- matrix(0, nrow = V, ncol = V)
    Omega.temp[lower.tri(Omega.temp, diag = T) & Omega == 1] <- c(1:((sum(Omega) - V)/2 + V))
    Omega.temp <- Omega.temp + t(Omega.temp) - diag(diag(Omega.temp))
    Q <- matrix(0, nrow = V^2, ncol = (sum(Omega)-V)/2 + V)
    Q[cbind(c(1:V^2), c(Omega.temp))] <- 1


    ## Setup Transpositon matrix such that vec(t(I-B)) = Tr %*% vec(I-B)
    Tr <- matrix(0, nrow = V^2, ncol = V ^2)
    Tr[matrix(c(1:V^2, (rep(c(1:V), V)-1) *V + rep(c(1:V), each = V)), ncol = 2)] <- 1

    ## precompute quantities of interest
    omega.inv <- solve(Omega.hat)
    eye.minus.b.inv <- solve(diag(rep(1, V)) - B.hat)
    omega.inv.eye.minus.b <- omega.inv %*% (diag(rep(1, V)) - B.hat)

    K <- J <-  matrix(0, nrow = sum(B) + (sum(Omega)-V)/2 + V,
                      ncol = sum(B) + (sum(Omega)-V)/2 + V)

    for(i in 1:dim(Y)[2]){
      l <- c(t(P) %*% c(omega.inv.eye.minus.b %*% (Y[,i] %*% t(Y[,i])) - t(eye.minus.b.inv)),
             t(Q) %*% c(omega.inv - omega.inv.eye.minus.b %*% Y[,i] %*% t(Y[,i]) %*% t(omega.inv.eye.minus.b)) / (-2) )
      K <- K + l %*% t(l)
    }
    K <- K / dim(Y)[2]

    # precompute term
    back_term <- omega.inv %*% (diag(rep(1, V)) - B.hat) %*% S %*% t((diag(rep(1, V)) - B.hat)) %*% omega.inv

    # information for beta elements
    J[1:sum(B), 1:sum(B)] <- t(P) %*% (S %x% omega.inv + (eye.minus.b.inv %x% t(eye.minus.b.inv)) %*% Tr) %*% P

    # information for omega elements (notice difference between this and expected info above)
    J[-c(1:sum(B)),-c(1:sum(B))] <- -1/2 * t(Q) %*%((omega.inv %x% omega.inv) - (omega.inv %x% back_term) -
                                                      (back_term %x% omega.inv)) %*% Q

    # co-information for omega elements
    J[c(1:sum(B)),-c(1:sum(B))] <- t(P) %*% ((S %*% t(diag(rep(1, V)) - B.hat) %*% omega.inv) %x% omega.inv) %*% Q
    J[-c(1:sum(B)),c(1:sum(B))] <- t(J[c(1:sum(B)),-c(1:sum(B))])

    # return estimate of variance
    return(solve(J) %*%K %*% t(solve(J)) / n)
  } else {
    stop("Invalid type given. Options are: observed, expected, sandwich")
  }

}


