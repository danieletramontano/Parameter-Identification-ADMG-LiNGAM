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
  weights_obs <- runif(p_edges, min = -w_lim, max = w_lim)
  # Generate random weights for hidden confounding
  weights_hid <- runif(2 * p_hid_edges, min = -w_lim, max = w_lim)
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

# Function that estimates weights for directed edges in a mixed graph for a LSEM
#
# Args:
#   data: A matrix where each row represents a data sample and each column represents a variable.
#   g_dir: A directed graph object (igraph).
#   g_bid: A bidirected graph object (igraph).
#   par_init: A numeric vector of initial values for the edge weights.
#   ker: A character string specifying the kernel to be used ("poly" for polynomial or "RBF" radial basis function).
#
# Returns:
#   A numeric vector of the optimized parameters.
parameter_estimation <- function(data,
                     g_dir,
                     g_bid,
                     par_init,
                     ker = "poly",
                     poly_degree = 2){
  dir_edges = matrix(as_edgelist(g_dir), ncol = 2)
  B = as.matrix(as_adjacency_matrix(g_dir))
  O = as.matrix(as_adjacency_matrix(g_bid))
  
  g_size = nrow(B)
  if(ker == "RBF"){
    adj_matrix <- matrix(0, nrow = g_size, ncol = g_size)
    for (i in 1:dim(dir_edges)[1]) {
      adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -par_init[i]
    }
    adj_matrix <- t(adj_matrix)
    
    # Calculate mix_matrix
    adj <- diag(g_size) + adj_matrix
    latent_data = t(adj %*%t(data))
    median_heu = apply(latent_data, 2, function (x) 1/sqrt(median(dist(x))/2))
  }
  fn  <- function(weight) {
    
    adj_matrix <- diag(g_size)
    for (i in 1:dim(dir_edges)[1]) {
      adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -weight[i]
    }
    adj_matrix <- t(adj_matrix)
    latent_data = t(adj_matrix %*%t(data))
    
    loss_fun = 0
    for(i in 1:(g_size-1)){
      for(j in (i+1):g_size){
        if(O[i, j] == 0){
          if(ker == "poly"){
            rbfkernel_1 <- polydot(degree = poly_degree)
            rbfkernel_2 <- polydot(degree = poly_degree)
          }else{
            rbfkernel_1 <- rbfdot(sigma = median_heu[i])
            rbfkernel_2 <- rbfdot(sigma = median_heu[j])
          }
          H1 = kernelMatrix(rbfkernel_1, matrix(c(latent_data[, i], latent_data[, i]), ncol = 2, byrow = F))
          H2 = kernelMatrix(rbfkernel_2, matrix(c(latent_data[, j], latent_data[, j]), ncol = 2, byrow = F))
          loss_fun = loss_fun + compute_HSIC(H1, H2)
        }
      }
    }
    if(ker == "poly"){
      loss_fun = log(loss_fun)
    }
    return(loss_fun)
  }
  
  op = optim(fn = fn, gr = function(weight) pracma::grad(fn, weight, heps=6e-06), par = par_init, method = "L-BFGS-B")
  return(op$par)
}



# Wrapper for the experiments
# Args:
#   df_params: A data frame containing parameters for repetitions, batch sizes,
#              polynomial degrees, and sample sizes.
#   g_bid: A bidirected graph object (igraph).
#   g_dir: A directed graph object (igraph).
#   par_bool: A logical value indicating whether to use parallel processing.
#   distribution: A character string specifying the distribution of the data
#                 (default: "Laplace").
#
# Returns:
#   A data frame containing the computed loss values for each method and parameter set.
#

fixed_graph_experiments <- function(df_params,
                                    g_bid,
                                    g_dir,
                                    par_bool = T,
                                    distribution = "Laplace"){
  
  B = as.matrix(as_adjacency_matrix(g_dir))
  O = as.matrix(as_adjacency_matrix(g_bid))
  g_size = nrow(B)
  if(dim(O)[1] == 0){
    O = matrix(0, nrow = g_size, ncol = g_size)
  }
  
  dir_edges = matrix(as_edgelist(g_dir), ncol = 2)
  n_dir_edges = length(dir_edges)/2
  
  g_size = nrow(B)
  bid_size = length(edges(g_bid))
  
  out_list = list()
  
  X = matrix(sample(c(-1, 1), 100*n_dir_edges, replace = T), ncol = n_dir_edges) * randomLHS(100, n_dir_edges)*3
  for (re in unique((df_params$reps))){
    out = NULL
    sampled_vars = mixed_graph_data_var(g_bid = g_bid,
                                        g_dir = g_dir)
    
    det_adj = 0
    while(abs(det_adj) < 1e-6){
      sampled_weights = mixed_graph_data_weights(g_bid = g_bid,
                                                 g_dir = g_dir)
      
      dir_edges = matrix(as_edgelist(g_dir), ncol = 2)
      
      edg_length <- dim(dir_edges)[1]
      adj <- diag(length(V(g_dir)))
      for (i in 1:edg_length) {
        adj[dir_edges[i, 1], dir_edges[i, 2]] <- -sampled_weights$weights_obs[i]
      }
      
      det_adj = det(adj)
    }
    sampled_data = mixed_graph_data(g_bid = g_bid,
                                    g_dir = g_dir, 
                                    data_size = max(df_params$n_sizes),
                                    distr = distribution,
                                    weights_obs = sampled_weights$weights_obs,
                                    weights_hid = sampled_weights$weights_hid,
                                    varr_obs = sampled_vars$varr_obs,
                                    varr_hid = sampled_vars$varr_hid)
    
    out$B = -sampled_data$adj
    out$Y = sampled_data$data
    out$weights = sampled_weights
    out$vars = sampled_vars
    
    out_list =  append(out_list, list(out))
  }
  
  if(par_bool == T){
    print("opened parallel threads")
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(
      n.cores,
      type = "FORK"
    )
    doParallel::registerDoParallel(cl = my.cluster)
  }
  
  (vals <- foreach(r = (df_params$reps),
                   poly_degree = (df_params$poly_degrees),
                   n = (df_params$n_sizes)) %dopar% {
                     
                     rep_ind = which(unique(df_params$reps) == r)
                     out = out_list[[rep_ind]]
                     true_weights = c()
                     
                     for(i in 1:dim(dir_edges)[1]){
                       true_weights[i] = out$B[dir_edges[i, 2], dir_edges[i, 1]]
                     }
                     
                     datas = out$Y
                     data = datas[1:n,]
                     
                     data =  data - t(matrix(colMeans(data), ncol(data), n))
                     data_cov = sample.cov(data)
                     
                     optim_fn <- function(par_init, 
                                          ker = "poly"){
                       if(ker == "RBF"){
                         adj_matrix <- matrix(0, nrow = g_size, ncol = g_size)
                         for (i in 1:dim(dir_edges)[1]) {
                           adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -par_init[i]
                         }
                         adj_matrix <- t(adj_matrix)
                         
                         # Calculate mix_matrix
                         adj <- diag(g_size) + adj_matrix
                         latent_data = t(adj %*%t(data))
                         median_heu = apply(latent_data, 2, function (x) 1/sqrt(median(dist(x))/2))
                       }
                       fn  <- function(weight) {
                         
                         adj_matrix <- diag(g_size)
                         for (i in 1:dim(dir_edges)[1]) {
                           adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -weight[i]
                         }
                         adj_matrix <- t(adj_matrix)
                         latent_data = t(adj_matrix %*%t(data))
                         
                         loss_fun = 0
                         for(i in 1:(g_size-1)){
                           for(j in (i+1):g_size){
                             if(O[i, j] == 0){
                               if(ker == "poly"){
                                 rbfkernel_1 <- polydot(degree = poly_degree)
                                 rbfkernel_2 <- polydot(degree = poly_degree)
                               }else{
                                 rbfkernel_1 <- rbfdot(sigma = median_heu[i])
                                 rbfkernel_2 <- rbfdot(sigma = median_heu[j])
                               }
                               H1 = kernelMatrix(rbfkernel_1, matrix(c(latent_data[, i], latent_data[, i]), ncol = 2, byrow = F))
                               H2 = kernelMatrix(rbfkernel_2, matrix(c(latent_data[, j], latent_data[, j]), ncol = 2, byrow = F))
                               loss_fun = loss_fun + compute_HSIC(H1, H2)
                             }
                           }
                         }
                         if(ker == "poly"){
                           loss_fun = log(loss_fun)
                         }
                         return(loss_fun)
                       }
                       
                       op = optim(fn = fn, gr = function(weight) pracma::grad(fn, weight, heps=6e-06), par = par_init, method = "L-BFGS-B")
                       return(op)
                     }
                     
                     #Regression
                     par_init = c()
                     for(i in 1:dim(dir_edges)[1]){
                       par_init[i] = data_cov[dir_edges[i, 1], dir_edges[i, 2]]/data_cov[dir_edges[i, 1], dir_edges[i, 1]]
                     }
                     
                     data_loss_cov = sum(((par_init - true_weights)^2/sum(true_weights^2)))
                     
                     #Optimize from regression coefficient
                     op = optim_fn(par_init, ker = "RBF")
                     data_loss_rbf <- sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     op = optim_fn(par_init, ker = "poly")
                     data_loss_poly <- sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     op = optim_fn(op$par, ker = "RBF")
                     data_loss_poly_rbf <- sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     #Empirical likelihood
                     
                     found_error = FALSE
                     # Use tryCatch to handle errors
                     found_error = tryCatch({
                       est <- sempl(Y = data, B = t(B), O = O)
                     }, error = function(e) {
                       # If an error occurs, set found_error to TRUE
                       found_error <- TRUE
                       return(found_error)
                     })
                     
                     if(length(found_error) == 1){
                       el_init = c(1:dim(dir_edges)[1])
                       data_loss_el =  sum(((el_init  - true_weights)^2/sum(true_weights^2)))
                     }else{
                       found_error = FALSE
                       el_init =  c()
                       for(i in 1:dim(dir_edges)[1]){
                         el_init[i] = est$B_hat[dir_edges[i, 2], dir_edges[i, 1]]
                       }
                       data_loss_el =  sum(((el_init  - true_weights)^2/sum(true_weights^2)))
                     }
                     
                     
                     #Optimize from empirical likelihood
                     op = optim_fn(el_init, ker = "RBF")
                     data_loss_el_optim_rbf = sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     op = optim_fn(el_init, ker = "poly")
                     data_loss_el_optim_poly = sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     op = optim_fn(el_init, ker = "RBF")
                     data_loss_el_optim_poly_rbf = sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     #Initialization at true value
                     op = optim_fn(true_weights, ker = "RBF")
                     data_loss_true_rbf =  sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     op = optim_fn(true_weights, ker = "poly")
                     data_loss_true_poly =  sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     op = optim_fn(op$par, ker = "RBF")
                     data_loss_true_poly_rbf =  sum(((op$par  - true_weights)^2/sum(true_weights^2)))
                     
                     res = list(reg_loss = data_loss_cov,
                                
                                reg_loss_optim_poly = data_loss_poly,
                                reg_loss_optim_poly_rbf = data_loss_poly_rbf,
                                reg_loss_optim_rbf = data_loss_rbf,
                                
                                tv_loss_optim_poly = data_loss_true_poly,
                                tv_loss_optim_rbf = data_loss_true_rbf,
                                tv_loss_optim_poly_rbf = data_loss_true_poly_rbf,
                                
                                el_loss = data_loss_el,
                                el_loss_optim_poly = data_loss_el_optim_poly,
                                el_loss_optim_rbf = data_loss_el_optim_rbf,
                                el_loss_optim_poly_rbf = data_loss_el_optim_poly_rbf,
                                found_error = found_error,
                                
                                repetition = r,
                                poly_degree = poly_degree,
                                sample_size = n)
                     return (res)
                   })
  
  df_reg_loss = df_reg_loss_optim_poly = df_reg_loss_optim_rbf = df_reg_loss_optim_poly_rbf = df_params
  df_tv_loss_optim_poly = df_tv_loss_optim_rbf = df_tv_loss_optim_poly_rbf = df_params
  df_el_loss = df_el_loss_optim_poly = df_el_loss_optim_rbf = df_el_loss_optim_poly_rbf= df_params
  
  df_reg_loss$method = 'REG'
  
  vals_length = length(vals)
  temp_vals = rep(0, vals_length)
  
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$reg_loss
  }
  df_reg_loss$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_reg_loss$found_error = temp_vals
  
  df_reg_loss_optim_poly$method = 'REG_poly_ker'
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$reg_loss_optim_poly
  }
  df_reg_loss_optim_poly$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_reg_loss_optim_poly$found_error = temp_vals
  
  df_reg_loss_optim_rbf$method = 'REG_rbf_ker'
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$reg_loss_optim_rbf
  }
  df_reg_loss_optim_rbf$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_reg_loss_optim_rbf$found_error = temp_vals
  
  df_reg_loss_optim_poly_rbf$method = 'REG_poly_rbf_ker'
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$reg_loss_optim_poly_rbf
  }
  df_reg_loss_optim_poly_rbf$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_reg_loss_optim_poly_rbf$found_error = temp_vals
  
  df_tv_loss_optim_poly$method = 'TV_poly_ker'
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$tv_loss_optim_poly
  }
  df_tv_loss_optim_poly$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_tv_loss_optim_poly$found_error = temp_vals
  
  df_tv_loss_optim_rbf$method = 'TV_rbf_ker'
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$tv_loss_optim_rbf
  }
  df_tv_loss_optim_rbf$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_tv_loss_optim_rbf$found_error = temp_vals
  
  df_tv_loss_optim_poly_rbf$method = 'TV_poly_rbf_ker'
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$tv_loss_optim_poly_rbf
  }
  df_tv_loss_optim_poly_rbf$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_tv_loss_optim_poly_rbf$found_error = temp_vals
  
  df_el_loss$method = 'EL'
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$el_loss
  }
  df_el_loss$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_el_loss$found_error = temp_vals
  
  df_el_loss_optim_poly$method = 'EL_poly_ker'
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$el_loss_optim_poly
  }
  df_el_loss_optim_poly$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_el_loss_optim_poly$found_error = temp_vals
  
  df_el_loss_optim_poly_rbf$method = 'EL_poly_rbf_ker'
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$el_loss_optim_poly_rbf
  }
  df_el_loss_optim_poly_rbf$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_el_loss_optim_poly_rbf$found_error = temp_vals
  
  df_el_loss_optim_rbf$method = 'EL_rbf_ker'
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$el_loss_optim_rbf
  }
  df_el_loss_optim_rbf$loss = temp_vals
  
  temp_vals = rep(0, vals_length)
  for(i in 1:vals_length){
    temp_vals[i] = vals[[i]]$found_error
  }
  df_el_loss_optim_rbf$found_error = temp_vals
  
  df <- rbind(df_reg_loss,
              df_reg_loss_optim_poly,
              df_reg_loss_optim_rbf,
              df_reg_loss_optim_poly_rbf,
              df_tv_loss_optim_poly,
              df_tv_loss_optim_rbf,
              df_tv_loss_optim_poly_rbf,
              df_el_loss,
              df_el_loss_optim_poly,
              df_el_loss_optim_rbf,
              df_el_loss_optim_poly_rbf)
  
  df$distribution = distribution
  
  if(par_bool == T){
    parallel::stopCluster(cl = my.cluster)
    print("closed parallel threads")
  }
  return(df)
}
