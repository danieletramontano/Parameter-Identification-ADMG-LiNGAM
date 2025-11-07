# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the current location of the file
source(file.path("setup.R"))



df_params <- expand.grid(
  n_sizes = c(10, 100, 300, 500, 700, 1000, 2000),
  poly_degrees = c(2),
  reps = c(1:10)
)


ww <-as.numeric(Sys.getenv("SLURM_PROCID"))+1; { ##Running on HPC using SLURM workload manager, comment if not using it
  #for(ww in c(1:5)){ ##Uncomment to run simulations in sequence
  if(ww == 1){
    set.seed(43)
    dc_sols = parallel_wrapper(df_params = df_params,
                                      g_bid = graph(edges = c(), directed = FALSE),
                                      g_dir = graph(edges = c(1, 2, 2, 3, 3, 1), directed = TRUE),
                                      distribution = "Laplace")

    dc_sols$graph = "3-Cycle"
    saveRDS(dc_sols,file="../Data/cycle_laplace_data_10_1.Rds")
  }
  if(ww == 2){
    set.seed(44)
    dc_sols = parallel_wrapper_el(df_params = df_params,
                                      g_bid = graph(edges = c(), directed = FALSE),
                                      g_dir = graph(edges = c(1, 2, 2, 3, 3, 1), directed = TRUE),
                                      distribution = "Laplace")

    dc_sols$graph = "3-Cycle"
    saveRDS(dc_sols,file="../Data/cycle_laplace_data_10_2.Rds")
  }

  if(ww == 3){
    set.seed(45)
    dc_sols = parallel_wrapper(df_params = df_params,
                                      g_bid = graph(edges = c(), directed = FALSE),
                                      g_dir = graph(edges = c(1, 2, 2, 3, 3, 1), directed = TRUE),
                                      distribution = "Laplace")

    dc_sols$graph = "3-Cycle"
    saveRDS(dc_sols,file="../Data/cycle_laplace_data_10_3.Rds")
  }
  if(ww == 4){
    set.seed(46)
    dc_sols = parallel_wrapper(df_params = df_params,
                                      g_bid = graph(edges = c(), directed = FALSE),
                                      g_dir = graph(edges = c(1, 2, 2, 3, 3, 1), directed = TRUE),
                                      distribution = "Laplace")

    dc_sols$graph = "3-Cycle"
    saveRDS(dc_sols,file="../Data/cycle_laplace_data_10_4.Rds")
  }

  if(ww == 5){
    set.seed(47)
    dc_sols = parallel_wrapper(df_params = df_params,
                                      g_bid = graph(edges = c(), directed = FALSE),
                                      g_dir = graph(edges = c(1, 2, 2, 3, 3, 1), directed = TRUE),
                                      distribution = "Laplace")

    dc_sols$graph = "3-Cycle"
    saveRDS(dc_sols,file="../Data/cycle_laplace_data_10_5.Rds")
  }
}