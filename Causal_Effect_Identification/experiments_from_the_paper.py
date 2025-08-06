#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 18:56:33 2024

@author: daniele
"""
import random
import identification



PARALLELIZE = True # set to False if don't want to parallelize
if PARALLELIZE:
    random.seed(10)
    results_n_25 = identification.parallel_rand_graphs_exp(n = 25)
    data_25 = identification.Countdata(n = 25,
                         reps = 5000,
                         counts = results_n_25,
                         ran = 100)
    identification.save_object(data_25)

    random.seed(10)
    results_n_50 = identification.parallel_rand_graphs_exp(n = 50)
    data_50 = identification.Countdata(n = 50,
                         reps = 5000,
                         counts = results_n_50,
                         ran = 100)
    identification.save_object(data_50)
else:
    random.seed(10)
    results_n_25 = identification.rand_graphs_exp(n = 25,
                                   ran = 50,
                                   rep = 50)
    data_25 = identification.Countdata(n = 25,
                         reps = 5000,
                         counts = results_n_25,
                         ran = 100)

    random.seed(10)
    results_n_50 = identification.rand_graphs_exp(n = 50)
    data_50 = identification.Countdata(n = 50,
                         reps = 5000,
                         counts = results_n_50,
                         ran = 10)


##If want to directly load the data uncomment the lines below
# data_25 = identification.load_object("Data/id_count-n_ 25 -medg_ 250-rep_5000.pickle")
# data_50 = identification.load_object("Data/id_count-n_ 50 -medg_ 500-rep_5000.pickle")

identification.prod_plot_mult_plots(files = [data_25, data_50],
                     label_fontsize=17,
                     subtitle_fontsize=25) # reproduce Fig. 10 of the paper