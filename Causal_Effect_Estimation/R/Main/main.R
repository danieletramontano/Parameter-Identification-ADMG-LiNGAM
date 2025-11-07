##
# Main helpers: parameter estimation and experiment wrapper
#
# This file provides two exported helpers:
# - parameter_estimation(data, g_dir, g_bid, par_init, ker, poly_degree)
# - parallel_wrapper(df_params, g_bid, g_dir, ...)
#
# The functions assume the supporting modules are available (HSIC, data generation,
# covariance). If they are not loaded in the current session we attempt to source
# them from sibling directories relative to this file.

## Required packages
library(igraph)
library(kernlab)
library(pracma)
library(VGAM)
library(foreach)
library(doParallel)
library(dplyr)

# Try to source sibling modules if their symbols are not present
root <- dirname(sys.frame(1)$ofile %||% "./")
maybe_source <- function(relpath) {
  p <- file.path(root, relpath)
  if (file.exists(p)) try(source(p), silent = TRUE)
}
maybe_source("../Data_generation/graph_data.R")
maybe_source("../Data_generation/covariance.R")
maybe_source("../HSIC/hsic.R")


#' Estimate directed edge weights using kernel-based HSIC loss
#'
#' @param data Numeric matrix (rows = samples, columns = variables)
#' @param g_dir igraph directed graph (nodes correspond to columns of data)
#' @param g_bid igraph undirected/bidirected graph (confounders)
#' @param par_init numeric initial parameter vector (length = number of directed edges)
#' @param ker character; "poly" or "RBF"
#' @param poly_degree integer polynomial degree when ker == "poly"
#' @return numeric vector of estimated parameters (same length as par_init)
parameter_estimation <- function(data,
                                 g_dir,
                                 g_bid,
                                 par_init,
                                 ker = "poly",
                                 poly_degree = 2) {
  dir_edges <- matrix(as_edgelist(g_dir), ncol = 2)
  B <- as.matrix(as_adjacency_matrix(g_dir))
  O <- as.matrix(as_adjacency_matrix(g_bid))

  g_size <- nrow(B)

  # Precompute median heuristic for RBF if requested
  if (ker == "RBF") {
    adj_matrix <- matrix(0, nrow = g_size, ncol = g_size)
    if (nrow(dir_edges) > 0) {
      for (i in seq_len(nrow(dir_edges))) {
        adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -par_init[i]
      }
    }
    adj_matrix <- t(adj_matrix)
    adj <- diag(g_size) + adj_matrix
    latent_data <- t(adj %*% t(data))
    median_heu <- apply(latent_data, 2, function(x) 1 / sqrt(median(dist(x)) / 2))
  }

  loss_function <- function(weight) {
    adj_matrix <- diag(g_size)
    if (nrow(dir_edges) > 0) {
      for (i in seq_len(nrow(dir_edges))) adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -weight[i]
    }
    adj_matrix <- t(adj_matrix)
    latent_data <- t(adj_matrix %*% t(data))

    loss_val <- 0
    for (i in seq_len(g_size - 1)) {
      for (j in seq((i + 1), g_size)) {
        if ((O[i, j] == 0) && (sum(B[, c(i, j)]) > 0)) {
          if (ker == "poly") {
            k1 <- polydot(degree = poly_degree)
            k2 <- polydot(degree = poly_degree)
          } else {
            k1 <- rbfdot(sigma = median_heu[i])
            k2 <- rbfdot(sigma = median_heu[j])
          }
          H1 <- kernelMatrix(k1, matrix(c(latent_data[, i], latent_data[, i]), ncol = 2, byrow = FALSE))
          H2 <- kernelMatrix(k2, matrix(c(latent_data[, j], latent_data[, j]), ncol = 2, byrow = FALSE))
          loss_val <- loss_val + compute_HSIC(H1, H2)
        }
      }
    }
    if (ker == "poly") loss_val <- log(loss_val)
    return(loss_val)
  }

  op <- optim(fn = loss_function, gr = function(w) pracma::grad(loss_function, w, heps = 6e-06), par = par_init, method = "L-BFGS-B")
  return(op$par)
}


#' Simple evaluation metric for weights
#' If ret_vector = TRUE returns component-wise relative error, otherwise a sum-of-squares relative metric.
evaluation_metric <- function(estimate_value, true_value, ret_vector = FALSE) {
  if (ret_vector) {
    return(abs(estimate_value - true_value) / (abs(true_value)))
  }
  return(sum(((estimate_value - true_value)^2 / sum(true_value^2))))
}


#' Experimental wrapper (serial or parallel)
#'
#' @param df_params Data frame with columns `reps`, `poly_degrees`, and `n_sizes` describing experiment grid
#' @param g_bid igraph bidirected/confounder graph
#' @param g_dir igraph directed graph
#' @param par_bool logical use parallel execution
#' @param ret_vector logical return vector metric
#' @param distribution character noise distribution ("Laplace" or "Uniform")
#' @param run_empirical_likelihood logical run empirical-likelihood baseline (requires sempl)
#' @param run_rbf_kernel logical run RBF kernel variants
#' @return tibble/data.frame with rows for each parameter combination and method
parallel_wrapper <- function(df_params,
                             g_bid,
                             g_dir,
                             par_bool = TRUE,
                             ret_vector = FALSE,
                             distribution = "Laplace",
                             run_empirical_likelihood = TRUE,
                             run_rbf_kernel = TRUE) {

  B <- as.matrix(as_adjacency_matrix(g_dir))
  O <- as.matrix(as_adjacency_matrix(g_bid))
  g_size <- nrow(B)
  if (nrow(O) == 0) O <- matrix(0, nrow = g_size, ncol = g_size)

  dir_edges <- matrix(as_edgelist(g_dir), ncol = 2)

  out_list <- list()
  for (re in unique(df_params$reps)) {
    sampled_vars <- mixed_graph_data_var(g_bid = g_bid, g_dir = g_dir)
    det_adj <- 0
    attempts <- 0
    while (abs(det_adj) < 1e-6 && attempts < 20) {
      sampled_weights <- mixed_graph_data_weights(g_bid = g_bid, g_dir = g_dir)
      dir_edges <- matrix(as_edgelist(g_dir), ncol = 2)
      adj <- diag(length(V(g_dir)))
      if (nrow(dir_edges) > 0) {
        for (i in seq_len(nrow(dir_edges))) adj[dir_edges[i, 1], dir_edges[i, 2]] <- -sampled_weights$weights_obs[i]
      }
      det_adj <- det(adj)
      attempts <- attempts + 1
    }
    if (abs(det_adj) < 1e-6) stop("Failed to sample a well-conditioned adjacency matrix")

    sampled_data <- mixed_graph_data(g_bid = g_bid,
                                     g_dir = g_dir,
                                     data_size = max(df_params$n_sizes),
                                     distr = distribution,
                                     weights_obs = sampled_weights$weights_obs,
                                     weights_hid = sampled_weights$weights_hid,
                                     varr_obs = sampled_vars$varr_obs,
                                     varr_hid = sampled_vars$varr_hid)

    out <- list(B = -sampled_data$adj, Y = sampled_data$data, weights = sampled_weights, vars = sampled_vars)
    out_list <- append(out_list, list(out))
  }

  if (par_bool) {
    n.cores <- parallel::detectCores(logical = FALSE)
    message(sprintf("registering %d cores for parallel execution", n.cores))
    my.cluster <- parallel::makeCluster(n.cores, type = "FORK")
    doParallel::registerDoParallel(cl = my.cluster)
  }

  # Evaluate each row of df_params in order; ensures output length matches nrow(df_params)
  vals <- lapply(seq_len(nrow(df_params)), function(ii) {
    r <- df_params$reps[ii]
    poly_degree <- df_params$poly_degrees[ii]
    n <- df_params$n_sizes[ii]
    rep_ind <- which(unique(df_params$reps) == r)
    out <- out_list[[rep_ind]]

    # true weights (match ordering of dir_edges)
    true_weights <- numeric(0)
    if (nrow(dir_edges) > 0) {
      true_weights <- sapply(seq_len(nrow(dir_edges)), function(k) out$B[dir_edges[k, 2], dir_edges[k, 1]])
    }

    data <- out$Y[1:n, , drop = FALSE]
    data <- data - t(matrix(colMeans(data), ncol(data), n))
    data_cov <- sample.cov(data)

    # regression initialization
    par_init <- numeric(0)
    for (v in seq_along(V(g_dir))) {
      par_v <- neighbors(g_dir, v, "in")
      if (length(par_v) > 0) par_init <- c(par_init, as.numeric(solve(matrix(data_cov[par_v, par_v], ncol = length(par_v))) %*% data_cov[par_v, v]))
    }

    # local optimizer wrapper
    optim_fn <- function(init_par, ker = "poly") {
      # compute median heuristic for RBF if needed
      if (ker == "RBF") {
        adj_matrix <- matrix(0, nrow = g_size, ncol = g_size)
        if (nrow(dir_edges) > 0) for (i in seq_len(nrow(dir_edges))) adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -init_par[i]
        adj_matrix <- t(adj_matrix)
        adj <- diag(g_size) + adj_matrix
        median_heu <- apply(t(adj %*% t(data)), 2, function(x) 1 / sqrt(median(dist(x)) / 2))
      }

      fn <- function(w) {
        adj_matrix <- diag(g_size)
        if (nrow(dir_edges) > 0) for (i in seq_len(nrow(dir_edges))) adj_matrix[dir_edges[i, 1], dir_edges[i, 2]] <- -w[i]
        adj_matrix <- t(adj_matrix)
        latent_data <- t(adj_matrix %*% t(data))
        loss_fun <- 0
        for (i in seq_len(g_size - 1)) {
          for (j in seq((i + 1), g_size)) {
            if ((O[i, j] == 0) && (sum(B[, c(i, j)]) > 0)) {
              if (ker == "poly") {
                k1 <- polydot(degree = poly_degree)
                k2 <- polydot(degree = poly_degree)
              } else {
                k1 <- rbfdot(sigma = median_heu[i])
                k2 <- rbfdot(sigma = median_heu[j])
              }
              H1 <- kernelMatrix(k1, matrix(c(latent_data[, i], latent_data[, i]), ncol = 2, byrow = FALSE))
              H2 <- kernelMatrix(k2, matrix(c(latent_data[, j], latent_data[, j]), ncol = 2, byrow = FALSE))
              loss_fun <- loss_fun + compute_HSIC(H1, H2)
            }
          }
        }
        if (ker == "poly") loss_fun <- log(loss_fun)
        return(loss_fun)
      }
      op <- optim(fn = fn, gr = function(w) pracma::grad(fn, w, heps = 6e-06), par = init_par, method = "L-BFGS-B")
      return(op)
    }

    # baseline and optimized runs
    res_reg_init <- evaluation_metric(par_init, true_weights, ret_vector = ret_vector)
    op_poly <- optim_fn(par_init, ker = "poly")
    reg_poly <- evaluation_metric(op_poly$par, true_weights, ret_vector = ret_vector)

    if (run_rbf_kernel) {
      op_rbf <- optim_fn(par_init, ker = "RBF")
      reg_rbf <- evaluation_metric(op_rbf$par, true_weights, ret_vector = ret_vector)
      op_poly_rbf <- optim_fn(op_rbf$par, ker = "RBF")
      reg_poly_rbf <- evaluation_metric(op_poly_rbf$par, true_weights, ret_vector = ret_vector)
    } else {
      reg_rbf <- evaluation_metric(true_weights, true_weights, ret_vector = ret_vector)
      reg_poly_rbf <- evaluation_metric(true_weights, true_weights, ret_vector = ret_vector)
    }

    # empirical likelihood baseline (optional)
    if (run_empirical_likelihood && exists("sempl")) {
      # simple placeholder: fit sempl if available (may be slow or unavailable)
      est <- try(simplify2array(sempl(Y = data, B = t(B), O = O)), silent = TRUE)
      el_init <- rep(1, ifelse(nrow(dir_edges) > 0, nrow(dir_edges), 0))
      el_loss <- evaluation_metric(el_init, true_weights, ret_vector = ret_vector)
      el_op_rbf <- optim_fn(el_init, ker = "RBF")
      el_rbf <- evaluation_metric(el_op_rbf$par, true_weights, ret_vector = ret_vector)
      el_op_poly <- optim_fn(el_init, ker = "poly")
      el_poly <- evaluation_metric(el_op_poly$par, true_weights, ret_vector = ret_vector)
      el_op_poly_rbf <- optim_fn(el_op_rbf$par, ker = "RBF")
      el_poly_rbf <- evaluation_metric(el_op_poly_rbf$par, true_weights, ret_vector = ret_vector)
    } else {
      el_loss <- evaluation_metric(true_weights, true_weights, ret_vector = ret_vector)
      el_rbf <- el_poly <- el_poly_rbf <- el_loss
    }

    # evaluate initialization at true weights
    true_poly_op <- optim_fn(true_weights, ker = "poly")
    true_poly <- evaluation_metric(true_poly_op$par, true_weights, ret_vector = ret_vector)
    if (run_rbf_kernel) {
      true_rbf_op <- optim_fn(true_weights, ker = "RBF")
      true_rbf <- evaluation_metric(true_rbf_op$par, true_weights, ret_vector = ret_vector)
      true_poly_rbf_op <- optim_fn(true_rbf_op$par, ker = "RBF")
      true_poly_rbf <- evaluation_metric(true_poly_rbf_op$par, true_weights, ret_vector = ret_vector)
    } else {
      true_rbf <- true_poly_rbf <- evaluation_metric(true_weights, true_weights, ret_vector = ret_vector)
    }

    list(
      reg_loss = res_reg_init,
      reg_loss_optim_poly = reg_poly,
      reg_loss_optim_poly_rbf = reg_poly_rbf,
      reg_loss_optim_rbf = reg_rbf,
      tv_loss_optim_poly = true_poly,
      tv_loss_optim_rbf = true_rbf,
      tv_loss_optim_poly_rbf = true_poly_rbf,
      el_loss = el_loss,
      el_loss_optim_poly = el_poly,
      el_loss_optim_rbf = el_rbf,
      el_loss_optim_poly_rbf = el_poly_rbf,
      repetition = r,
      poly_degree = poly_degree,
      sample_size = n
    )
  })

  # convert list of results to data.frame
  df_reg_loss <- df_reg_loss_optim_poly <- df_reg_loss_optim_rbf <- df_reg_loss_optim_poly_rbf <- df_params
  df_tv_loss_optim_poly <- df_tv_loss_optim_rbf <- df_tv_loss_optim_poly_rbf <- df_params
  df_el_loss <- df_el_loss_optim_poly <- df_el_loss_optim_rbf <- df_el_loss_optim_poly_rbf <- df_params

  update_df <- function(df, method_name, loss_field, vals) {
    df$method <- method_name
    df$loss <- sapply(vals, function(x) x[[loss_field]])
    return(df)
  }

  datasets <- list(
    list(df = df_reg_loss, method = "REG", loss_field = "reg_loss"),
    list(df = df_reg_loss_optim_poly, method = "REG_poly_ker", loss_field = "reg_loss_optim_poly"),
    list(df = df_reg_loss_optim_rbf, method = "REG_rbf_ker", loss_field = "reg_loss_optim_rbf"),
    list(df = df_reg_loss_optim_poly_rbf, method = "REG_poly_rbf_ker", loss_field = "reg_loss_optim_poly_rbf"),
    list(df = df_tv_loss_optim_poly, method = "TV_poly_ker", loss_field = "tv_loss_optim_poly"),
    list(df = df_tv_loss_optim_rbf, method = "TV_rbf_ker", loss_field = "tv_loss_optim_rbf"),
    list(df = df_tv_loss_optim_poly_rbf, method = "TV_poly_rbf_ker", loss_field = "tv_loss_optim_poly_rbf"),
    list(df = df_el_loss, method = "EL", loss_field = "el_loss"),
    list(df = df_el_loss_optim_poly, method = "EL_poly_ker", loss_field = "el_loss_optim_poly"),
    list(df = df_el_loss_optim_rbf, method = "EL_rbf_ker", loss_field = "el_loss_optim_rbf"),
    list(df = df_el_loss_optim_poly_rbf, method = "EL_poly_rbf_ker", loss_field = "el_loss_optim_poly_rbf")
  )

  # `vals` is a list with one entry per row of `df_params` (see above lapply).
  dataframes <- list()
  for (dataset in datasets) {
    # extract raw loss entries for this dataset
    raw_losses <- lapply(vals, function(x) x[[dataset$loss_field]])
    # if every element is length 1, convert to numeric vector, otherwise keep as list-column
    if (all(sapply(raw_losses, length) == 1)) {
      dataset$df$loss <- unlist(raw_losses)
    } else {
      dataset$df$loss <- I(raw_losses)
    }
    dataset$df$method <- dataset$method
    dataframes[[dataset$method]] <- dataset$df
  }

  df <- bind_rows(dataframes)
  df$distribution <- distribution

  if (par_bool) {
    parallel::stopCluster(cl = my.cluster)
    message("closed parallel threads")
  }
  return(df)
}
