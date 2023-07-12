#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import calculate_polarization, stars_polarization
import os
import matplotlib.pyplot as plt
import numpy as np

star_name = 'HD14069'

base_path = os.path.join('..', '..', 'Low polarized stars', star_name, 'reduced', star_name )
csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
pol_dict = calculate_polarization(csv_file_name)

for filter in pol_dict.keys():
    plt.plot(pol_dict[filter][0], pol_dict[filter][1], 'o-', label=f'filter {filter}')

plt.axhline(stars_polarization[star_name], color='r', linestyle='--', label='Literature')
plt.xlabel('Time (MJD)')
plt.ylabel('Polarization (%)')
plt.title(f'Polazation as a function of time for the star {star_name}')
plt.legend()
plt.show()

