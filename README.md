# Code for Reproducing the Experiments in "PARAMETER IDENTIFICATION IN LINEAR NON-GAUSSIAN CAUSAL MODELS UNDER GENERAL CONFOUNDING"

## Overview

This repository contains the code for reproducing the experiments presented in the paper "PARAMETER IDENTIFICATION IN LINEAR NON-GAUSSIAN CAUSAL MODELS UNDER GENERAL CONFOUNDING".

The code is organized into three main parts:

1. **Causal Effect Identification (Python)**

2. **Causal Effect Estimation (R)**

3. **Model Equivalence (Macaulay2)**

## 1. Causal Effect Identification (Python)

This folder contains the Python code and necessary resources to certify the criterion of Theorems 4.2 -4.3, and to reproduce the experiments of Section 8.1.

### Contents:

- **requirements.yml**: YML file for installing the required Python libraries.
- **example_notebook_identification.ipynb**: Jupyter notebook demonstrating how to use the algorithms on a simple example.
- **identification.py**: Python script that contains the functions to implement the algorithms.
- **experiment_code.py**: Python script to fully reproduce the experiments from the paper.
- **Data/**: Directory containing the output data of the experiments.

### Instructions:

1. **Setup Environment**:
   - Install Anaconda or Miniconda.

   - Navigate to the folder and create the environment using the YML file:

     ```bash
     conda env create -f environment.yml
     conda activate causal_effect_identification
     ```

2. **Run Example Notebook**:
   - Open `example_notebook.ipynb` in Jupyter Notebook or JupyterLab to see a step-by-step guide on using the algorithms.

3. **Reproduce Experiments**:
   - Execute `experiment_code.py` to reproduce the experiments and generate output data.

## 2. Causal Effect Estimation (R)

This folder contains the R code and necessary resources to execute the optimization described in Lemma 8.1, and to reproduce the experiments of Section 8.2.

### Contents:

- **install_packages.R**: R script for installing the required R packages.
- **example_notebook.Rmd**: R Markdown notebook demonstrating how to use the algorithms on a simple example.
- **experiment_code.R**: R script to fully reproduce the experiments from the paper.
- **output_data/**: Directory containing the output data of the experiments.

### Instructions:

1. **Setup Environment**:
   - Ensure you have R and RStudio installed.

   - Install the required packages by running:

     ```r
     source("install_packages.R")
     ```

2. **Run Example Notebook**:
   - Open `example_notebook.Rmd` in RStudio to see a step-by-step guide on using the algorithms.

3. **Reproduce Experiments**:
   - Execute `iv_simulations.R` to reproduce the experiments and generate output data for the IV graph.
   - Execute `dc_simulations.R` to reproduce the experiments and generate output data for the graph in Fig. 11.
   - Execute `cycle_simulations_laplace.R` to reproduce the experiments and generate output data for the *3-cycle* of Fig. 9, and the error terms samples from a Laplace distribution.
   - Execute `cycle_simulations_uniform.R` to reproduce the experiments and generate output data for the *3-cycle* of Fig. 9, and the error terms samples from an Uniform distribution.
   - Execute `plots.R` to reproduce Fig. 12 from the main paper, and Fig.1-2 from the supplement.

## 3. Model Equivalence (Macaulay2)

This folder contains the Macaulay2 code for the model identification algorithm described in Theorem 6.2, and to reproduce the experiments from Example 6.1.

### Contents:

- **model_identification.m2**: Macaulay2 script containing the code for the algorithm described in Theorem 6.2.
- **example_script.m2**: Macaulay2 script to with a simple example of how to use the algorithm.
- **experiment_from_the_paper.m2**: Macaulay2 script to replicate the reproduce the experiments from Example 6.1.

### Instructions:

1. **Setup Environment**:
   - Ensure you have Macaulay2 installed. Download it from [Macaulay2 official website](http://www2.macaulay2.com/Macaulay2/).

2. **Run the Script**:
   - Open Macaulay2 and load the script:

     ```m2
     load "model_identification.m2"
     ```

   - Follow the instructions within the script to replicate the experiments.