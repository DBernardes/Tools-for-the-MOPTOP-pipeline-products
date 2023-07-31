#!/usr/bin/env python

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np
from tools import track_obj_over_images


star_name = 'BD+32 3739'
min, max =  59775, 60051
experiment = 'several positions in image'
src_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star_name, experiment, 'cumulative', star_name)
fig, axs = plt.subplots(2, 3, figsize=(15, 7), sharex=True, sharey=True)
axs[0, 0].set_ylabel('Y axis - camera 1')
axs[1, 0].set_ylabel('Y axis - camera 2')

for camera in [3, 4]:
    for idx, filter in enumerate(['V']):
        ax = axs[camera-3, idx]
        new_path = os.path.join(src_path, f'cam{camera}')
        xcoord, ycoord, _ = track_obj_over_images(new_path)
        ax.plot(xcoord, ycoord, 'o-', alpha = 0.2)

        medx, medy = np.median(xcoord), np.median(ycoord)
        # ax.plot(medx, medy, '*k')
        # ax.annotate(f'({medx:.0f}, {medy:.0f})', (medx*1.1,medy*1.1), fontsize=10, ha='left')
        
        ax.axvline(512, color='r', linestyle='--', alpha=0.75)
        ax.axhline(512, color='r', linestyle='--', alpha=0.75)
        if camera == 3:
            ax.set_title(f'Filter {filter}')
        else:
            ax.set_xlabel('X axis')


plt.xlim(0, 1024)
plt.ylim(0, 1024)
plt.show()