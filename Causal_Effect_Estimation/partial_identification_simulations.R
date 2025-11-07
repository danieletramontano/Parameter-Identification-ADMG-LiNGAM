# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the current location of the file
source(file.path("packages_install.R"))
source(file.path("utils.R"))


df_params <- expand.grid(
  n_sizes = c(100, 500),
  poly_degrees = c(2),
  reps = c(1:5)
)


refactor_df <- function(df, edge_list =  c("2->4", "3->4")) {
  loss_list <- df$loss
  num_dims <- length(loss_list[[1]])

  # Initialize an empty list to collect new data frames
  df_list <- list()
  # Iterate over dimensions to process the loss
  for (i in 1:num_dims) {
    temp_df <- df %>%
      mutate(
        loss = sapply(loss_list, function(x) x[i]),
        graph = edge_list[i]
      )

    # Add the processed data frame to the list
    df_list[[i]] <- temp_df
  }

  # Combine all the data frames
  final_df <- bind_rows(df_list)
  return(final_df)
}


ww <-as.numeric(Sys.getenv("SLURM_PROCID"))+1; { ##Running on HPC using SLURM workload manager, comment if not using it
  #for(ww in c(1:2)){ ##Uncomment to run simulations in sequence
  if(ww == 1){
    set.seed(42)
    pi_sols = fixed_graph_experiments(df_params = df_params,
                                      g_bid = graph(edges = c(1, 2,
                                                              2, 4,
                                                              3,4), directed = FALSE),
                                      g_dir = graph(edges = c(2, 4,
                                                              3, 4), directed = TRUE),
                                      ret_vector = T,
                                      distribution = "Laplace",
                                      run_empirical_likelyhood = FALSE,
                                      run_rbf_kernel = FALSE)

    pi_sols_refactor = refactor_df(pi_sols)
    saveRDS(pi_sols_refactor,file="Data/pi_laplace_data.Rds")
  }
  if(ww == 2){
    set.seed(42)
    pi_sols = fixed_graph_experiments(df_params = df_params,
                                      g_bid = graph(edges = c(1, 2,
                                                              1, 4,
                                                              2, 4,
                                                              2, 3), directed = FALSE),
                                      g_dir = graph(edges = c(1, 4,
                                                              2, 4,
                                                              3, 4), directed = TRUE),
                                      ret_vector = T,
                                      distribution = "Uniform",
                                      run_empirical_likelyhood = FALSE,
                                      run_rbf_kernel = FALSE)

    pi_sols_refactor = refactor_df(pi_sols)
    saveRDS(pi_sols_refactor,file="Data/pi_uniform_data.Rds")
  }
}