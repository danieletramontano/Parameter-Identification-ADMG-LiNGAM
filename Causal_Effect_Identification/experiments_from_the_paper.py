#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 18:56:33 2024

@author: daniele
"""

from identification import parallel_ranggraph, rangraphsimul, count_data, prod_plot_mult_plots, save_object, load_object
import random



run_par = True # set to False if don't want to parallelize
if run_par is True:
    random.seed(10)
    results_n_25 = parallel_ranggraph(n = 25)
    data_25 = count_data(n = 25,
                         reps = 5000,
                         counts = results_n_25,
                         ran = 10)
    save_object(data_25)
    
    
    random.seed(10)
    results_n_50 = parallel_ranggraph(n = 50)
    data_50 = count_data(n = 50,
                         reps = 5000,
                         counts = results_n_50,
                         ran = 10)
    save_object(data_50)
else:
    random.seed(10)
    results_n_25 = rangraphsimul(n = 25)
    data_25 = count_data(n = 25,
                         reps = 5000,
                         counts = results_n_25,
                         ran = 10)
    
    random.seed(10)
    results_n_50 = rangraphsimul(n = 50)
    data_50 = count_data(n = 50,
                         reps = 5000,
                         counts = results_n_50,
                         ran = 10)


##If want to directly load the data uncomment the line below
# data_25 = load_object("Data/id_count-n_ 25 -medg_ 250-rep_5000.pickle")
# data_50 = load_object("Data/id_count-n_ 50 -medg_ 500-rep_5000.pickle")

prod_plot_mult_plots(files = [data_25, data_50],
                     label_fontsize=17,
                     subtitle_fontsize=25,
                     legend_fontsize=18) # reproduce Fig. 10
