#!/usr/bin/env python

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import track_obj_over_images


star_name = 'BD+32 3739'
src_path = os.path.join('..', 'Low polarized stars', star_name, star_name, 'selected_files', 'cam3')

sep = 0
for filter in ['V', 'R', 'I']:
    new_path = os.path.join(src_path, f'{filter} Filter')
    xcoord, ycoord = track_obj_over_images(new_path)
    plt.plot(xcoord, ycoord-sep, 'o-', alpha = 0.2, label=f'{filter}')
    sep+=150



plt.xlim(0, 1024)
plt.ylim(0, 1024)
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('sensor of the CCD')
plt.legend()
plt.show()