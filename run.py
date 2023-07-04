#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import manipulate_csv_file, calculate_polarization
import os
import matplotlib.pyplot as plt
import numpy as np


base_path = os.path.join('..', 'Low polarized stars', 'HD14069', 'reduced', 'HD14069' )
csv_file_name = os.path.join(base_path, 'raw_data.csv')
manipulate_csv_file(csv_file_name)
csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
pol_dict = calculate_polarization(csv_file_name)

for filter in pol_dict.keys():
    # plt.errorbar(pol_dict[filter]['mjd'], pol_dict[filter]['pol'], pol_dict[filter]['std'], fmt='o-', label=f'filter {filter}')
    pol = np.asarray(pol_dict[filter]['pol'])*100
    plt.plot(pol_dict[filter]['mjd'], pol, 'o-', label=f'filter {filter}')

plt.axhline(0.111, color='r', linestyle='--', label='Literature')
plt.xlabel('Time (MJD)')
plt.ylabel('Polarization (%)')
plt.title('Polazation as a function of time for the star HD14069')
plt.legend()
plt.show()