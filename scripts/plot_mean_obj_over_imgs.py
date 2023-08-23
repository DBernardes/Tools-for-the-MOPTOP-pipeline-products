#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np
from tools import (
    track_obj_over_images,
    get_obj_coords,
    calculate_mean_box_pixels,
    sort_files,
)
from astropy.coordinates import SkyCoord
import astropy.units as u
import astropy.io.fits as fits
from sys import exit

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
# RA, DEC = "19:03:33.2", "40:53:16.5" ORIGINAL
# RA, DEC = "19:03:34.9", "40:53:07.1" INTERESTING
# RA, DEC = "19:03:33.24", "40:53:21.5"
RA, DEC = "19:03:36.26", "40:53:26.53"


# fig, axs = plt.subplots(1, 1, figsize=(15, 7), sharex=True, sharey=True)
for camera in [3, 4]:
    images_list = sort_files(src_path, f"{camera}_e")
    mean_pixels, mjds = [], []
    for image in images_list:
        file_name = os.path.join(src_path, image)
        x, y, mjd = get_obj_coords(file_name, RA, DEC)
        x = int(x)
        y = int(y)
        tmp = calculate_mean_box_pixels(file_name, x, y, 50)
        mean_pixels.append(tmp)
        mjds.append((mjd))
    plt.plot(mjds, mean_pixels, "o", alpha=0.5, label=f"camera {camera}")
plt.xlabel("MJD")
plt.ylabel("Photons")
plt.legend()
plt.show()
