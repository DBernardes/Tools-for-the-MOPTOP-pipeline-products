#!/usr/bin/env python

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import track_obj_over_images
import numpy as np
import pandas as pd

min, max =  59768, 59825
star_name = 'BD+32 3739'
base_path = os.path.join('..', '..', 'Low polarized stars', star_name)
fig, axs = plt.subplots(3, 3, figsize=(15, 7))
camera = 4

for idx, filter in enumerate(['V', 'R', 'I']):
    new_path  = os.path.join(base_path, star_name, 'selected_files', f'{min}-{max}', f'cam{camera}', f'{filter} Filter')
    xcoord, ycoord, mjd = track_obj_over_images(new_path)
    ax = axs[0, idx]
    ax.plot(mjd, xcoord, 'o', label=f'x coord')
    ax.plot(mjd, ycoord, 'o', label=f'y coord')
    ax.set_ylabel('Pixels')
    ax.set_xlim(mjd[0], mjd[-1])
    ax.legend()
    ax.set_title(f'Filter {filter}')

    new_path  = os.path.join(base_path, 'reduced', star_name, 'manipulated_data.csv')
    df = pd.read_csv(new_path)
    rows = df.loc[df['wave'] == f'MOP-{filter}']
    photons = rows[f's1_src_cam{camera-2}']
    mean = np.median(photons)
    std = np.std(np.abs(photons - mean))
    ax = axs[1, idx]
    ax.plot(rows['mjd'], photons, 'o', label=f'cam{camera}')
    ax.set_ylabel('Photons')
    ax.set_xlim(mjd[0], mjd[-1])
    ax.set_ylim(mean-1.5*std, mean+1.5*std)
    ax.invert_yaxis()
    ax.legend()

    ax = axs[2, idx]
    q = rows['q_avg']
    u = rows['u_avg']
    tmp = np.mean(q-u)
    ax.plot(rows['mjd'], q, 'o', label=f'q')
    ax.plot(rows['mjd'], u+tmp, 'o', label=f'u + {tmp:.2f}')
    ax.set_xlabel('Time (MJD)')
    ax.set_xlim(mjd[0], mjd[-1])
    ax.legend()


plt.show()