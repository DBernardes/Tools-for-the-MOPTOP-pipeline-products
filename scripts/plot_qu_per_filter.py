#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import sort_qu_per_filter, low_polarized_stars, high_polarized_stars
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sys import exit


fig, axs = plt.subplots(2, 5, figsize=(18, 10), sharex=True, sharey='row')
for ax in axs[1]: 
    ax.set_xlabel('Time (MJD)')

i=0
for stokes_param in ['q', 'u']:
    axs[i,0].set_ylabel(f'{stokes_param} value')
    for idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
        ax = axs[i, idx]
        if i == 0:
            ax.set_title(f'Filter {filter}')
        colors = ['b', 'r', 'g', 'k', 'm', 'c', 'y']
        sep = 0
        for key in high_polarized_stars.keys():
            color = colors.pop(0)
            base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'High polarized stars', key, 'reduced', key )
            csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
            qu_dict = sort_qu_per_filter(csv_file_name)
            mjd = qu_dict[filter]['mjd']
            qu_val = np.asarray(qu_dict[filter][stokes_param]) - sep
            qu_err = np.asarray(qu_dict[filter][f'{stokes_param}_err'])
            ax.errorbar(mjd, qu_val, yerr=qu_err, fmt=f'{color}o', alpha=0.25, label=f'{key}')
            #ax.set_ylim(-0.1,  0.1)
            #ax.set_xlim(59600, 60100)
            ax.legend()
            sep += 0.01
    i+=1
plt.show()