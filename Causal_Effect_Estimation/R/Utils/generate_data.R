# Calculate sample covariance matrix
#
# This function calculates the sample covariance matrix of two datasets.
# If only one dataset is provided, it calculates the sample covariance matrix of that dataset.
# If two datasets are provided, it calculates the covariance matrix between them.
#
# Args:
#   X: A numeric matrix or data frame.
#   Y: An optional numeric matrix or data frame with the same number of rows as X.
#
# Returns:
#   A covariance matrix.
#
sample.cov <- function(X, Y=NULL){
  n <- nrow(X)
  Xr <- X - t(matrix(colMeans(X), ncol(X), n))
  if (is.null(Y)){
    return(t(Xr)%*%Xr/n)
  } else {
    if (nrow(Y) != n) stop("X and Y must have the same number of rows.")
    Yr <- Y - t(matrix(colMeans(Y), ncol(Y), n))
    return(t(Xr)%*%Yr/n)
  }
}

runiff <- function(n, min = 0, max = 1){
  if(sample(c(1,-1), 1) > 0){
    return(runif(n, min = 0.5, max = max))
  }
  return(runif(n, min = min, max = -0.5))
}

# Generate random weights for mixed graph data
#
# This function generates random weights for the adjacency matrix of a mixed graph.
# The weights are generated separately for directed and bidirected graphs.
#
# Args:
#   g_bid: A bidirected graph object (igraph).
#   g_dir: A directed graph object (igraph).
#   w_lim: A numeric value specifying the range of weights.
#
# Returns:
#   A list containing weights for observed edges and hidden confounding.
#
mixed_graph_data_weights <- function(g_bid = NULL, g_dir = NULL, w_lim = 5) {
  # Count the edges in the directed and bidirected graphs
  p_edges <- dim(matrix(as_edgelist(g_dir), ncol = 2))[1]
  p_hid_edges <- dim(matrix(as_edgelist(g_bid), ncol = 2))[1]
  # Generate random weights for observed edges
  weights_obs <- runiff(p_edges, min = -w_lim, max = w_lim)
  # Generate random weights for hidden confounding
  weights_hid <- runiff(2 * p_hid_edges, min = -w_lim, max = w_lim)
  # Return the results
  return(list(weights_obs = weights_obs, weights_hid = weights_hid))
}

# Generate random variance for mixed graph data
#
# This function generates random variances for the exogenous noises in a mixed graph.
# Variances can be the same for all variables or different, and can be scaled by a specified ratio.
#
# Args:
#   g_bid: A bidirected graph object (igraph).
#   g_dir: A directed graph object (igraph).
#   same_var: A logical value indicating whether to use the same variance for all variables.
#   stn_ratio: A numeric value specifying the signal-to-noise ratio for hidden confounding.
#
# Returns:
#   A list containing variances for observed variables and hidden confounding.
#
mixed_graph_data_var <- function(g_bid = NULL, g_dir = NULL, same_var = FALSE, stn_ratio = 1) {
  # Generate random variance for the exogenous noises.

  g_size <- length(V(g_dir))

  p_hid_edges <- length(as_edgelist(g_bid))

  if (same_var) {
    varr_obs <- rep(1, g_size)
    varr_hid <- rep(1, p_hid_edges)
  } else {
    varr_obs <- runif(g_size, 0.2, 3)
    varr_hid <- rep(0, p_hid_edges)
    if(p_hid_edges > 0){
      for (i in 1:p_hid_edges) {
        varr_hid[i] <- varr_obs[(as_edgelist(g_bid))[i][[1]]] * stn_ratio
      }
    }
  }

  return(list(varr_obs = varr_obs, varr_hid = varr_hid))
}


# Generate data based on mixed graph structure
#
# This function generates data based on the structure of a mixed graph, including both directed and bidirected edges.
# The data can be generated from various distributions and includes options for specifying weights and variances.
#
# Args:
#   g_bid: A bidirected graph object (igraph).
#   g_dir: A directed graph object (igraph).
#   data_size: An integer specifying the number of data points to generate.
#   distr: A character string specifying the distribution of the data (default: 'Laplace').
#   weights_obs: A numeric vector of weights for observed edges.
#   weights_hid: A numeric vector of weights for hidden confounding.
#   varr_obs: A numeric vector of variances for observed variables.
#   varr_hid: A numeric vector of variances for hidden confounding.
#   conf_type: An integer specifying the type of confounding (default: 2).
#
# Returns:
#   A list containing generated data, adjacency matrix, data without confounding, and latent errors.
#
mixed_graph_data <- function(g_bid = NULL, g_dir = NULL, data_size = 500, distr = 'Laplace',
                             weights_obs = NULL, weights_hid = NULL, varr_obs = NULL, varr_hid = NULL,
                             conf_type = 2) {
  dir_edges = matrix(as_edgelist(g_dir), ncol = 2)

  edg_length <- dim(dir_edges)[1]

  # Create the adjacency matrix
  adj <- diag(length(V(g_dir)))
  for (i in 1:edg_length) {
    adj[dir_edges[i, 1], dir_edges[i, 2]] <- -weights_obs[i]
  }
  adj <- t(adj)

  # Calculate B matrix
  mix_matrix <- solve(adj)
  err <- matrix(0, nrow = data_size, ncol = length(V(g_dir)))

  for (j in 1:length(V(g_dir))) {
    if (distr == 'Laplace') {
      err[, j] <- rlaplace(data_size, 0, 1) * varr_obs[j]
    } else if(distr == 'Uniform'){
      err[, j] <- runif(data_size, -1, 1) * varr_obs[j]
    }
  }

  err_no_conf <- err
  data_no_conf <- t(mix_matrix %*% t(err_no_conf))
  data_no_conf <- data_no_conf - t(matrix(colMeans(data_no_conf), ncol(data_no_conf), data_size))

  bid_edges <- matrix(as_edgelist(g_bid), ncol = 2)

  if(dim(bid_edges)[1]>0){
    for(reps in 1:2){
      for (j in 1:dim(bid_edges)[1]) {
        conf_err <- numeric(data_size)
        #conf_type <- sample(0:1, 1)
        for (i in 1:data_size) {
          if (distr == 'Laplace') {
            conf_err[i] <- rlaplace(1, 0, 1) * varr_hid[j]
          } else if(distr == 'Uniform'){
            conf_err[i] <- runif(1, -1, 1) * varr_hid[j]
          }
        }

        if (conf_type == 1) {
          err[, bid_edges[j, 1]] <- err[, bid_edges[j, 1]] * conf_err * weights_hid[2 * j - 1]
          err[, bid_edges[j, 2]] <- err[, bid_edges[j, 2]] + weights_hid[2 * j] * conf_err^2
        } else {
          err[, bid_edges[j, 1]] <- err[, bid_edges[j, 1]] + weights_hid[2 * j -1 ] * conf_err
          err[, bid_edges[j, 2]] <- err[, bid_edges[j, 2]] + weights_hid[2 * j] * conf_err
        }
      }
    }
  }

  # Generate data
  data <- t(mix_matrix %*% t(err))
  data <- data - t(matrix(colMeans(data), ncol(data), data_size))


  return(list(data = data, adj = adj, data_no_conf = data_no_conf, lat_err = err))
}
