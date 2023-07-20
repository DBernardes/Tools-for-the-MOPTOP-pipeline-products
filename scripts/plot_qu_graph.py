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
from circle_fit import taubinSVD

# rows, columns = 2, 3
# fig, axs = plt.subplots(rows, columns, figsize=(20, 4), sharex=True, sharey=True)
# for idxy in range(rows):
#     axs[idxy,0].set_ylabel('u')
#     for idxx in range(columns):
#         ax = axs[idxy, idxx]
#         star = list(low_polarized_stars.keys())[idxx+columns*idxy]
#         colors = ['b', 'r', 'g', 'y', 'm']
#         base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star, 'reduced', star )
#         csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
#         qu_dict = sort_qu_per_filter(csv_file_name)
#         for idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
#             color = colors.pop(0)
#             q = qu_dict[filter]['q']
#             u = qu_dict[filter]['u']
#             q, u = sigma_clipping(q, u, 5, 1)
#             meanq, meanu = np.mean(q), np.mean(u)
        
#             ax.plot(q, u, f'{color}o', alpha=0.25, label=f'{filter}')
#             ax.axhline(0, color='r', linestyle='--', alpha=0.25)
#             ax.axvline(0, color='r', linestyle='--', alpha=0.25)
#             ax.plot(meanq, meanu, f'k*')
#             ax.annotate(f'({meanq:.3f},{meanu:.3f})', (meanq*0.95,meanu), fontsize=10, ha='right')
#             if idxy == rows-1:
#                 ax.set_xlabel('q')
#             ax.legend()
#             ax.grid()
#             ax.set_title(star)

# plt.show()


#--------------------------------------------------------------------------------------------------------------------------

base_path = os.path.join('..', '..', 'Pol charact MOPTOP')
csv_file_name = os.path.join(base_path, 'Low polarized stars', 'mean_qu_values.csv')
df = pd.read_csv(csv_file_name)
for star in df['star']:
    if star not in ['GD319', 'HD14069', 'BD+32 3739']:
        df.drop(df[df['star'] == star].index, inplace = True)


stars = high_polarized_stars.keys()
len_stars = len(stars)
fig, axs = plt.subplots(len_stars, 5, sharey='row', sharex='row', figsize=(5*len_stars, 10))
for row_idx, star in enumerate(stars):
    row_axs = axs[row_idx]
    colors = ['b', 'r', 'g', 'k', 'm']
    csv_file_name = os.path.join(base_path, 'High polarized stars', star, 'reduced', star , 'manipulated_data.csv')
    qu_dict = sort_qu_per_filter(csv_file_name)
    for col_idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
        ax = row_axs[col_idx]
        color = colors.pop(0)
        q, u = qu_dict[filter]['q'], qu_dict[filter]['u']
        rows = df.loc[df['filter'] == filter]
        inst_q, inst_u = np.mean(rows['q']), np.mean(rows['u'])
        q-=inst_q
        u-=inst_u

        point_coordinates = [(q_val, u_val) for (q_val, u_val) in zip(q, u)]
        xc, yc, r, sigma = taubinSVD(point_coordinates)
        circle1 = plt.Circle((xc, yc), r, color=color, fill=False)
        ax.add_patch(circle1)
        ax.plot(xc, yc, f'{color}+')

        ax.plot(q, u, f'{color}o', alpha=0.25, label=f'({xc:.3f},{yc:.3f},{r:.3f})')
        ax.axhline(0, color='r', linestyle='--', alpha=0.25)
        ax.axvline(0, color='r', linestyle='--', alpha=0.25)
        ax.legend(fontsize=8, loc='upper right')
        ax.grid()
        if row_idx == 0:
            ax.set_title(f'filter {filter}')
        if row_idx == len_stars-1:
            ax.set_xlabel(f'q')
        if col_idx == 0:
            ax.set_ylabel(f'Star {star}\nu')
        
plt.show()