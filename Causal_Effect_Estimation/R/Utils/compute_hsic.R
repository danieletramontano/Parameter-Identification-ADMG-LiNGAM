# Compute Hilbert-Schmidt Independence Criterion (HSIC)
#
# This function computes the Hilbert-Schmidt Independence Criterion (HSIC) between two kernel matrices.
#
# Args:
#   K: A square kernel matrix.
#   L: A square kernel matrix of the same dimensions as K.
#
# Returns:
#   A numeric value representing the HSIC.
#
compute_HSIC <- function(K, L){
  dim_k = dim(K)

  if (length(dim_k)!= 2 | dim_k[1] != dim_k[2]){
    stop("K must be a square matrix")
  }
  else dim_k = dim_k[1]

  dim_l = dim(L)
  if (length(dim_l)!= 2 | dim_l[1] != dim_l[2]){
    stop("L must be a square matrix")
  }
  else dim_l = dim_l[1]
  if (dim_l != dim_k){
    stop("K and L must have the same dimension")
  }
  else n = dim_l

  H = diag(n)-1/n
  M = K %*% H %*% L %*% H

  h = sum(diag(M))/((n-1)**2)
  return(h)
}