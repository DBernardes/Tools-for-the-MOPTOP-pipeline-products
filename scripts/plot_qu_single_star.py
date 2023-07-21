#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import sort_qu_per_filter
import os
import matplotlib.pyplot as plt
import numpy as np
from sys import exit


star_name = 'BD+32 3739'
fig, axs = plt.subplots(5, 1, figsize=(18, 5), sharex=True)
axs[0].set_title(f'Star {star_name}')

for idx, filter in enumerate(['B', 'V', 'R', 'I', 'L']):
    ax = axs[idx]
    #ax.set_title(f'Filter {filter}')
    base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star_name, 'reduced', star_name )
    csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
    qu_dict = sort_qu_per_filter(csv_file_name)
    q = np.asarray(qu_dict[filter]['q'])
    u = np.asarray(qu_dict[filter]['u'])
    ax.plot(qu_dict[filter]['mjd'], q, 'bo', alpha=0.5, label=f'q')
    ax.plot(qu_dict[filter]['mjd'], u, 'ro', alpha=0.5, label=f'u')
    ax.set_ylim(-0.1,  0.1)
    #ax.set_xlim(59600, 60100)
    ax.legend()
    ax.set_ylabel(f'Filter {filter}')
plt.xlabel('Time (MJD)')
plt.show()