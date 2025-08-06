# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the current location of the file
source(file.path("packages_install.R"))
source(file.path("utils.R"))


df_params <- expand.grid(
  n_sizes = c(10, 100, 300, 500, 700, 1000, 2000),
  poly_degrees = c(2,3,4,5,6),
  reps = c(1:50)
)

#ww <-as.numeric(Sys.getenv("SLURM_PROCID"))+1; { ##Running on HPC using SLURM workload manager, comment if not using it

set.seed(42)
dc_sols = fixed_graph_experiments(df_params = df_params,
                                  g_bid = make_graph(edges = c(1, 2, 1,3), directed = FALSE),
                                  g_dir = make_graph(edges = c(1,2, 1, 3, 2, 3), directed = TRUE),
                                  distribution = "Laplace",
                                  run_empirical_likelyhood = FALSE,
                                  only_polynomial_kernel = TRUE)

dc_sols$graph = 'Double Confounder'
saveRDS(dc_sols,file="Data/dc_laplace_data_hyperp.Rds")
