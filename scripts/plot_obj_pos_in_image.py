#!/usr/bin/env python

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import track_obj_over_images


star_name = 'BD+32 3739'
min, max =  59768, 59825
src_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star_name, 'several positions in image', star_name)


new_path = os.path.join(src_path, f'{filter} Filter')
xcoord, ycoord, _ = track_obj_over_images(src_path, '3_e')
plt.plot(xcoord, ycoord, 'ro-', alpha = 0.2, label='cam3')

xcoord, ycoord, _ = track_obj_over_images(src_path, '4_e')
plt.plot(xcoord, ycoord, 'bo-', alpha = 0.2, label='cam4')


plt.hlines(512, 0, 1024, color='r', linestyle='--', alpha=0.25)
plt.vlines(512, 0, 1024, color='r', linestyle='--', alpha=0.25)
plt.xlim(0, 1024)
plt.ylim(0, 1024)
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.title('CCD sensor')
plt.legend()
plt.show()