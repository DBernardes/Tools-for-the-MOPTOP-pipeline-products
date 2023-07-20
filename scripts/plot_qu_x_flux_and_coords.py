#!/usr/bin/env python

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import get_coords_in_series
import numpy as np
import pandas as pd
from scipy import stats

min, max =  59860,59932
star_name = 'VICyg12'
camera = 4
alpha = 0.7
fontsize = 9
base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'High polarized stars', star_name)
fig, axs = plt.subplots(3, 3, figsize=(15, 7), sharey='row')
axs[2,0].set_xlabel(f'Photons cam{camera}')
axs[2,1].set_xlabel('X axis')
axs[2,2].set_xlabel('Y axis')


def plot_graph(ax, x, q, u, res_q, res_u, filter, tmp):
    ax.plot(x, q, 'o', label=f'q', alpha=alpha)
    ax.plot(x, u, 'o', label=f'u+{tmp:.2f}', alpha=alpha)
    ax.annotate(f'coef_q={res_q.statistic:.2f}, pval_q={res_q.pvalue:.2f}\ncoef_u={res_u.statistic:.2f}, pval_u={res_u.pvalue:.2f}',
            xy=(.05, .95), xycoords='axes fraction',
            ha='left', va='top',
            fontsize=fontsize,
            bbox=dict(boxstyle="round",
                      fc="lightgray", ec="black", alpha=0.5))
    ax.set_ylabel(f'{filter} filter')
    ax.legend(loc='lower right')
    return

def calc_plot_parameters(df, x_str):
    x = df[x_str]
    q = df['q_avg']
    u = df['u_avg']
    res_q = stats.spearmanr(x, q)
    res_u = stats.spearmanr(x, u)
    return x, q, u, res_q, res_u

for idx, filter in enumerate(['V','R', 'I']):
    print(f'processing the images of filter {filter}...')

    new_path  = os.path.join(base_path, 'reduced', star_name, 'manipulated_data.csv')
    df = pd.read_csv(new_path)
    rows = df.loc[df['wave'] == f'MOP-{filter}']
    rows = rows.loc[min <= rows['mjd']]
    rows = rows.loc[rows['mjd'] <= max]
    tmp= np.mean(rows['q_avg']-rows['u_avg'])
    rows['u_avg'] += tmp
    
    rows = rows.sort_values(by=[f's1_src_cam{camera-2}'])
    ax = axs[idx, 0]
    photons, q, u, res_p1q, res_p1u = calc_plot_parameters(rows, f's1_src_cam{camera-2}')
    plot_graph(ax, photons, q, u, res_p1q, res_p1u, filter, tmp)

    #------------------------------------------------------------
    new_path  = os.path.join(base_path, star_name, 'selected_files', f'{min}-{max}', f'cam{camera}', f'{filter} Filter')
    xcoord, ycoord = get_coords_in_series(new_path, rows['date'], rows['mjd'])
    rows['xcoord'] = xcoord
    rows['ycoord'] = ycoord

    ax = axs[idx, 1]
    rows = rows.sort_values(by=[f'xcoord'])
    coord, q, u, res_xq, res_xu = calc_plot_parameters(rows, f'xcoord')
    plot_graph(ax, coord, q, u, res_xq, res_xu, filter, tmp)

    ax = axs[idx, 2]
    rows = rows.sort_values(by=[f'ycoord'])
    coord, q, u, res_yq, res_yu = calc_plot_parameters(rows, f'ycoord')
    plot_graph(ax, coord, q, u, res_yq, res_yu, filter, tmp)
    
plt.show()


