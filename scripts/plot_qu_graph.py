#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import sort_qu_per_filter, low_polarized_stars, high_polarized_stars, sigma_clipping
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sys import exit

rows, columns = 2, 3
fig, axs = plt.subplots(rows, columns, figsize=(20, 4), sharex=True, sharey=True)
for idxy in range(rows):
    axs[idxy,0].set_ylabel('u')
    for idxx in range(columns):
        ax = axs[idxy, idxx]
        star = list(low_polarized_stars.keys())[idxx+columns*idxy]
        colors = ['b', 'r', 'g', 'y', 'm']
        base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star, 'reduced', star )
        csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
        qu_dict = sort_qu_per_filter(csv_file_name)
        for idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
            color = colors.pop(0)
            q = qu_dict[filter]['q']
            u = qu_dict[filter]['u']
            q, u = sigma_clipping(q, u, 5, 1)
            meanq, meanu = np.mean(q), np.mean(u)
            ax.plot(q, u, f'{color}o', alpha=0.25, label=f'{filter}')
            ax.axhline(0, color='r', linestyle='--', alpha=0.25)
            ax.axvline(0, color='r', linestyle='--', alpha=0.25)
            ax.plot(meanq, meanu, f'k*')
            ax.annotate(f'({meanq:.3f},{meanu:.3f})', (meanq*0.95,meanu), fontsize=10, ha='right')
            if idxy == rows-1:
                ax.set_xlabel('q')
            ax.legend()
            ax.grid()
            ax.set_title(star)
plt.show()


#--------------------------------------------------------------------------------------------------------------------------

# rows, columns = 2,3
# fig, axs = plt.subplots(rows, columns, figsize=(20, 4), sharex='col', sharey='row')
# for idxy in range(rows):
#     for idxx in range(columns):
#         ax = axs[idxy, idxx]
#         ax.set_ylabel('u')
#         try:
#             star = list(high_polarized_stars.keys())[idxx+columns*idxy]
#         except:
#             break
#         colors = ['b', 'r', 'g', 'k', 'm']
#         base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'High polarized stars', star, 'reduced', star )
#         csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
#         qu_dict = sort_qu_per_filter(csv_file_name)
#         for idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
#             color = colors.pop(0)
#             ax.plot(qu_dict[filter]['q'], qu_dict[filter]['u'], f'{color}o', alpha=0.25, label=f'{filter}')
#             ax.axhline(0, color='r', linestyle='--', alpha=0.25)
#             ax.axvline(0, color='r', linestyle='--', alpha=0.25)
#             if idxy == rows-1:
#                 ax.set_xlabel('q')
#             ax.legend()
#             ax.grid()
#             ax.set_title(star)
# plt.show()