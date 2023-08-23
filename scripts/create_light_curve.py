#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np
from tools import (
    get_obj_coords,
    sort_files,
)
import astropy.io.fits as fits

star_name = "GRB 230818A"
experiment = "all data"
src_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
    experiment,
    star_name,
)

RA, DEC = "19:03:40", "+40:52:35"
size = 30
light_curve = []
image_files = os.listdir(src_path)
image_files = sort_files(src_path, f"{3}_e")

for file in image_files:
    print(file)
    file = os.path.join(src_path, file)
    x, y, _ = get_obj_coords(file, RA, DEC)
    max_val = calculate_box_pixels(file, x, y, size=size, calc_method="max")
    x, y = give_val_get_xy_coords(file, x, y, max_val, size=size)
    try:
        fwhm = calc_FWHM(file, x, y)
        photons = calc_psf_photons(file, x, y, 3 * fwhm)
        light_curve.append(photons)
    except Exception:
        continue

plt.plot(light_curve, "ob")
plt.show()
