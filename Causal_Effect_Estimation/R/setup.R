setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to the current location of the file
is_installed <- function(package) {
  return(package %in% installed.packages()[,"Package"])
}

# Install missing packages
install_missing_packages <- function(packages) {
  for (pkg in packages) {
    if (!is_installed(pkg)) {
      install.packages(pkg)
    }
  }
}



if(!is_installed("empLikSem")){
  install.packages("empLikSem_1.0.1.tar", repos=NULL, type='source')
}
library(empLikSem)


# List of required packages
packages <- c('pracma',
              'distr',
              'LaplacesDemon',
              'lhs',
              'kernlab',
              'foreach',
              'pcalg',
              'igraph',
              'doParallel',
              'ggplot2',
              'scales',
              'tidyverse')

# Function to check if package is installed


# Install missing packages
install_missing_packages(packages)

# Load required libraries
library(pracma)
library(distr)
library(LaplacesDemon)
library(lhs)
library(kernlab)
library(foreach)
library(pcalg)
library(igraph)
library(doParallel)
library(tidyverse)
library(dplyr)
library(tidyr)

# Now, the required libraries are installed and loaded

source(file.path("Main/main.R"))
