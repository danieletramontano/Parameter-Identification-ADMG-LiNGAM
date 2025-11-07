# Generates a plot of mean loss against sample size for different methods and graphs.
# Arguments:
#   data: (data frame) Input data containing loss values and parameters.
#   poly_degree: (numeric) Polynomial degree for loss calculation. Default is 2.
#   dist_type: (character) Type of noise distribution. Default is "Laplace".
#   inits: (character vector) List of initialization methods. Default is c("EL", "REG", "TV").
#   ker_type: (character) Type of kernel. Default is "poly".
#   show_graphs: (character vector) Graphs to include in the plot. Default is c('Instrumental Variable', 'Double Confounder', '3-Cycle').
#   title_size: (numeric) Font size for plot title. Default is 15.
#   axis_size: (numeric) Font size for axis labels. Default is 13.
#   legend_size: (numeric) Font size for legend text. Default is 15.
#   grid_size: (numeric) Font size for grid text. Default is 15.
#   figure_name: (character) Name of the saved figure file. Default is "" (empty string).
#   save_fig: (logical) Indicates whether to save the figure. Default is FALSE.

# Returns:
#   ggplot object representing the plot.
produce_plot = function(data,
                        poly_degree = 2,
                        dist_type = "Uniform",
                        inits = list("EL", "REG", "TV"),
                        ker_type = "poly",
                        show_graphs = c('Instrumental Variable', 'Double Confounder', '3-Cycle'),
                        title_size = 15,
                        axis_size = 13,
                        legend_size = 15,
                        grid_size = 15,
                        figure_name = "",
                        save_fig = F,
                        loss = loss,
                        log_trans = T){

  linetypes_methods = c(
    "REG_poly_ker"= "longdash",
    "REG_poly_rbf_ker"= "longdash",
    "REG_rbf_ker"= "longdash",

    "TV_poly_ker"  = "dashed",
    "TV_poly_rbf_ker"  = "dashed",
    "TV_rbf_ker"  = "dashed",

    "EL_poly_ker"  = "dotted",
    "EL_poly_rbf_ker"  = "dotted",
    "EL_rbf_ker"  = "dotted",

    "EL" = "solid"
  )


  methods = lapply(inits, FUN = function(x) paste0(x,"_",ker_type,"_ker" ))
  #methods = append(methods, "EL")

  data$graph_o = factor(data$graph, levels=show_graphs)

  #data = filter(data, method %in% methods, poly_degrees == poly_degree, graph %in% show_graphs)
  data = filter(data, method %in% methods, graph %in% show_graphs)
  if(!is.null(dist_type)){
    data = filter(data, distribution == dist_type)
  }

  data$loss = as.numeric(data$loss)

  fig = ggplot(data,aes(poly_degrees,
                        y = loss,
                        linetype = factor(method))) +
    stat_summary(fun = mean,
                 geom = "line") +
    scale_linetype_manual(values = linetypes_methods) +
    xlab("Sample Size") +
    ylab("Mean Loss") +
    # facet_grid(~poly_degrees) +
    labs(title = paste(dist_type, "Noise")) +
    theme_light() +
    theme(plot.title = element_text(face = "bold", size = title_size),
          legend.title = element_blank(),
          legend.position = "bottom",
          axis.text = element_text(size = axis_size),
          text = element_text(size = legend_size),
          strip.text.x = element_text(size = grid_size),
          aspect.ratio = 1,
          axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = grid::unit(c(0, 4, 0, 2), "mm"))
  if (log_trans) {
    fig <- fig + scale_y_continuous(trans = "log10")
  }
  if(save_fig == T){
    ggsave(plot = fig, dpi = 300, filename = paste0(getwd(),"/Figures/", figure_name))
  }
  fig
  return(fig)
}
produce_boxplot = function(data,
                        poly_degree = 2,
                        dist_type = "Uniform",
                        inits = list("EL", "REG", "TV"),
                        ker_type = "poly",
                        show_graphs = c('Instrumental Variable', 'Double Confounder', '3-Cycle'),
                        title_size = 15,
                        axis_size = 13,
                        legend_size = 15,
                        grid_size = 15,
                        figure_name = "",
                        save_fig = F,
                        loss = loss,
                        log_trans = T){

  linetypes_methods = c(
    "REG_poly_ker"= "longdash",
    "REG_poly_rbf_ker"= "longdash",
    "REG_rbf_ker"= "longdash",

    "TV_poly_ker"  = "dashed",
    "TV_poly_rbf_ker"  = "dashed",
    "TV_rbf_ker"  = "dashed",

    "EL_poly_ker"  = "dotted",
    "EL_poly_rbf_ker"  = "dotted",
    "EL_rbf_ker"  = "dotted",

    "EL" = "solid"
  )


  methods = lapply(inits, FUN = function(x) paste0(x,"_",ker_type,"_ker" ))

  data$graph_o = factor(data$graph, levels=show_graphs)

  data = filter(data, method %in% methods, graph %in% show_graphs)
  if(!is.null(dist_type)){
    data = filter(data, distribution == dist_type)
  }

  data$loss = as.numeric(data$loss)

    fig = ggplot(data, aes(x = factor(poly_degrees),
                         y = loss,
                         fill = factor(method))) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    scale_fill_grey(start = 0.3, end = 0.8) +
    facet_grid(~graph) +
    xlab("Degree of the Polynomial Kernel") +
    ylab("Loss") +
    labs(title = paste(dist_type, "Noise")) +
    theme_light() +
    theme(plot.title = element_text(face = "bold", size = title_size),
          legend.title = element_blank(),
          legend.position = "bottom",
          axis.text = element_text(size = axis_size),
          text = element_text(size = legend_size),
          strip.text.x = element_text(size = grid_size),
          aspect.ratio = 1,
          axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = grid::unit(c(0, 4, 0, 2), "mm"))
  if (log_trans) {
    fig <- fig + scale_y_continuous(trans = "log10")
  }
  if(save_fig == T){
    ggsave(plot = fig, dpi = 300, filename = paste0(getwd(),"/Figures/", figure_name))
  }
  fig
  return(fig)
}


# Merges multiple datasets stored as RDS files in a directory into a single dataset.
# Arguments:
#   data_directory: (character) Path to the directory containing the RDS files.
#                   If NULL, the function looks for files in the "Data/" directory of the current working directory.
#   save_data:      (logical) Indicates whether to save the merged dataset as an RDS file. Default is FALSE.
#   data_name:      (character) Name of the merged dataset file. Default is "complete_data.Rds".
# Returns:
#   Merged dataset as a data frame.
merge_datasets <- function(data_directory = NULL,
                           save_data = F,
                           data_name = "complete_data.Rds"){
  if(is.null(directory_path)){
    directory_path <- paste0(getwd(),"/Data/")
  }

  file_list <- list.files(path = directory_path, full.names = TRUE)
  data_list <- c()
  for(i in c(1:length(file_list))){
    data_list[[i]] <- readRDS(file_list[[i]])
  }

  merged_data <- do.call(rbind, data_list)
  if(save_data){
    saveRDS(merged_data,file=paste0(directory_path,data_name))
  }
  return(merged_data)
}

### If want to merge the datasets in a folder uncomment the following line
# complete_data =  merge_datasets(data_directory = NULL,
#                                 save_data = F,
#                                 data_name = "complete_data.Rds")


complete_data = readRDS("Data/complete_data.Rds")
produce_plot(complete_data,
             poly_degree = 2,
             dist_type = "Laplace",
             inits = list("EL", "REG", "TV"),
             ker_type = "poly",
             title_size = 27,
             axis_size = 23.5,
             legend_size = 27,
             grid_size = 27,
             figure_name = "laplace_poly.png",
             save_fig = F)

produce_plot(complete_data,
             poly_degree = 2,
             dist_type = "Uniform",
             inits = list("EL", "REG", "TV"),
             ker_type = "poly",
             title_size = 27,
             axis_size = 23.5,
             legend_size = 27,
             grid_size = 27,
             figure_name = "uniform_poly.png",
             save_fig = F)

produce_plot(complete_data,
             poly_degree = 2,
             dist_type = "Laplace",
             inits = list("EL", "REG", "TV"),
             ker_type = "rbf",
             title_size = 27,
             axis_size = 23.5,
             legend_size = 27,
             grid_size = 27,
             figure_name = "laplace_rbf.png",
             save_fig = F)

produce_plot(complete_data,
             poly_degree = 2,
             dist_type = "Uniform",
             inits = list("EL", "REG", "TV"),
             ker_type = "rbf",
             title_size = 27,
             axis_size = 23.5,
             legend_size = 27,
             grid_size = 27,
             figure_name = "unform_rbf.png",
             save_fig = F)

produce_plot(complete_data,
             poly_degree = 2,
             dist_type = "Uniform",
             inits = list("REG", "TV"),
             show_graphs = c('2 -> 4', '3 -> 4'),
             ker_type = "rbf",
             title_size = 27,
             axis_size = 23.5,
             legend_size = 27,
             grid_size = 27,
             figure_name = "pi_uniform.png",
             save_fig = F)

produce_plot(complete_data,
             poly_degree = 2,
             dist_type = "Laplace",
             inits = list("REG", "TV"),
             show_graphs = c('2 -> 4', '3 -> 4'),
             ker_type = "rbf",
             title_size = 27,
             axis_size = 23.5,
             legend_size = 27,
             grid_size = 27,
             figure_name = "pi_uniform.png",
             save_fig = F)

dc_hyperparameters_data = readRDS("Data/dc_hyperparameters.Rds")
produce_boxplot(dc_hyperparameters_data,
             poly_degree = 2,
             dist_type = "Laplace",
             inits = list("REG", "TV"),
             show_graphs = c('Double Confounder'),
             ker_type = "rbf",
             title_size = 27,
             axis_size = 23.5,
             legend_size = 27,
             grid_size = 27,
             figure_name = "pi_uniform.png",
             save_fig = F)