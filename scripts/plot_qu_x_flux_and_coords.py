#!/usr/bin/env python

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import get_coords_in_series
import numpy as np
import pandas as pd
from scipy import stats

min, max =  0,1e10
star_name = 'BD+32 3739'
experiment = 'several positions in image'
camera = 4
alpha = 0.7
fontsize = 9
base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star_name, experiment, 'cumulative')
fig, axs = plt.subplots(3, 3, figsize=(15, 7), sharey='row')
axs[2,0].set_xlabel(f'Photons cam{camera}')
axs[2,1].set_xlabel('X axis')
axs[2,2].set_xlabel('Y axis')


def plot_graph(ax, x, q, q_err, u, u_err, res_q, res_u, filter, tmp):
    ax.errorbar(x, q, yerr=q_err, fmt='o', label=f'q', alpha=alpha)
    ax.errorbar(x, u, yerr=u_err, fmt='o', label=f'u+{tmp:.2f}', alpha=alpha)
    ax.annotate(f'coef_q={res_q.statistic:.3f}, pval_q={res_q.pvalue:.3f}\ncoef_u={res_u.statistic:.3f}, pval_u={res_u.pvalue:.3f}',
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
    q_err = df['q_err']
    u_err = df['u_err']
    res_q = stats.spearmanr(x, q)
    res_u = stats.spearmanr(x, u)
    return x, q, q_err, u, u_err, res_q, res_u

for idx, filter in enumerate(['V']):
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
    photons, q, q_err, u, u_err, res_p1q, res_p1u = calc_plot_parameters(rows, f's1_src_cam{camera-2}')
    plot_graph(ax, photons, q, q_err, u, u_err, res_p1q, res_p1u, filter, tmp)

    #------------------------------------------------------------
    new_path  = os.path.join(base_path, star_name, f'cam{camera}')#, 'selected_files', f'{min}-{max}', f'cam{camera}', f'{filter} Filter')
    xcoord, ycoord = get_coords_in_series(new_path, rows['date'], rows['mjd'], f'{camera}_e')
    rows['xcoord'] = xcoord
    rows['ycoord'] = ycoord

    ax = axs[idx, 1]
    rows = rows.sort_values(by=[f'xcoord'])
    coord, q, q_err, u, u_err, res_xq, res_xu = calc_plot_parameters(rows, f'xcoord')
    plot_graph(ax, coord, q, q_err, u, u_err, res_xq, res_xu, filter, tmp)

    ax = axs[idx, 2]
    rows = rows.sort_values(by=[f'ycoord'])
    coord, q, q_err, u, u_err, res_yq, res_yu = calc_plot_parameters(rows, f'ycoord')
    plot_graph(ax, coord, q, q_err, u, u_err, res_yq, res_yu, filter, tmp)
    
plt.show()


