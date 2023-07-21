#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import sort_qu_per_filter, low_polarized_stars, high_polarized_stars, sigma_clipping, read_calculated_qu_values, LIMIT_MJD
import os
import matplotlib.pyplot as plt
import pandas as pd
from circle_fit import taubinSVD

calculated_qu = read_calculated_qu_values()
base_path = os.path.join('..', '..', 'Pol charact MOPTOP')

rows, columns = 2, 3
colors = ['b', 'r', 'g', 'y', 'm']
fig, axs = plt.subplots(rows, columns, figsize=(20, 4), sharex=True, sharey=True)

for idxy in range(rows):
    axs[idxy,0].set_ylabel('u')
    for idxx in range(columns):
        ax = axs[idxy, idxx]
        star = list(low_polarized_stars.keys())[idxx+columns*idxy]
        csv_file_name = os.path.join(base_path, 'Low polarized stars', star, 'reduced', star, 'manipulated_data.csv')

        df = pd.read_csv(csv_file_name)
        df = df.loc[df['mjd'] > LIMIT_MJD]
        for idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
            color = colors[idx]
            df_rows = df.loc[df['wave'] == f'MOP-{filter}']
            q = df_rows['q_avg']
            u = df_rows['u_avg']
            q, u = sigma_clipping(q, u)
            inst_q, inst_u = calculated_qu[filter]
        
            ax.plot(q, u, f'{color}o', alpha=0.25, label=f'{filter}')
            ax.axhline(0, color='r', linestyle='--', alpha=0.25)
            ax.axvline(0, color='r', linestyle='--', alpha=0.25)
            ax.plot(inst_q, inst_u, f'k*')
            ax.annotate(f'({inst_q:.3f},{inst_u:.3f})', (inst_q*0.95,inst_u), fontsize=10, ha='right')
            if idxy == rows-1:
                ax.set_xlabel('q')
            ax.legend()
            ax.grid()
            ax.set_title(star)

plt.show()


#--------------------------------------------------------------------------------------------------------------------------




# stars = high_polarized_stars.keys()
# len_stars = len(stars)
# fig, axs = plt.subplots(len_stars, 5, sharey='row', sharex='row', figsize=(5*len_stars, 10))
# for row_idx, star in enumerate(stars):
#     row_axs = axs[row_idx]
#     colors = ['b', 'r', 'g', 'k', 'm']
#     csv_file_name = os.path.join(base_path, 'High polarized stars', star, 'reduced', star , 'manipulated_data.csv')
#     qu_dict = sort_qu_per_filter(csv_file_name)
#     for col_idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
#         ax = row_axs[col_idx]
#         color = colors.pop(0)
#         q, u = qu_dict[filter]['q'], qu_dict[filter]['u']
#         inst_q, inst_u = calculated_qu[filter]
#         q-=inst_q
#         u-=inst_u

#         point_coordinates = [(q_val, u_val) for (q_val, u_val) in zip(q, u)]
#         xc, yc, r, sigma = taubinSVD(point_coordinates)
#         circle1 = plt.Circle((xc, yc), r, color=color, fill=False)
#         ax.add_patch(circle1)
#         ax.plot(xc, yc, f'{color}+')

#         ax.plot(q, u, f'{color}o', alpha=0.25, label=f'({xc:.3f},{yc:.3f},{r:.3f})')
#         ax.axhline(0, color='r', linestyle='--', alpha=0.25)
#         ax.axvline(0, color='r', linestyle='--', alpha=0.25)
#         ax.legend(fontsize=8, loc='upper right')
#         ax.grid()
#         if row_idx == 0:
#             ax.set_title(f'filter {filter}')
#         if row_idx == len_stars-1:
#             ax.set_xlabel(f'q')
#         if col_idx == 0:
#             ax.set_ylabel(f'Star {star}\nu')
        
# plt.show()