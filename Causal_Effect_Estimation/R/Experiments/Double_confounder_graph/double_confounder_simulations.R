# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the current location of the file
source(file.path("setup.R"))



df_params <- expand.grid(
  n_sizes = c(10, 100, 300, 500, 700, 1000, 2000),
  poly_degrees = c(2),
  reps = c(1:50)
)


ww <-as.numeric(Sys.getenv("SLURM_PROCID"))+1; { ##Running on HPC using SLURM workload manager, comment if not using it
  #for(ww in c(1:2)){ ##Uncomment to run simulations in sequence
  if(ww == 1){
    set.seed(42)
    dc_sols = parallel_wrapper(df_params = df_params,
                                      g_bid = graph(edges = c(1, 2, 1,3), directed = FALSE),
                                      g_dir = graph(edges = c(1,2, 1, 3, 2, 3), directed = TRUE),
                                      distribution = "Laplace")

    dc_sols$graph = 'Double Confounder'
    saveRDS(dc_sols,file="../Data/dc_laplace_data.Rds")
  }
  if(ww == 2){
    set.seed(42)
    dc_sols = parallel_wrapper(df_params = df_params,
                                      g_bid = graph(edges = c(1, 2, 1,3), directed = FALSE),
                                      g_dir = graph(edges = c(1,2, 1, 3, 2, 3), directed = TRUE),
                                      distribution = "Uniform")

    dc_sols$graph = 'Double Confounder'
    saveRDS(dc_sols,file="../Data/dc_uniform_data.Rds")
  }
}
