#!/usr/bin/env python


__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
import pandas as pd
from tools import get_obj_coords
import astropy.units as u
from matplotlib import transforms

hsize = 1024
RA, DEC = "19:03:33.2", "40:53:16.5"
star_name = "GRB 230818A"
experiment = "all data"
base_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
)

fits_file = os.path.join(base_path, "all data", star_name, "3_e_20230818_6_3_4_1.fits")
data = fits.getdata(fits_file)
obj_x, obj_y, mjd = get_obj_coords(fits_file, RA, DEC)
x, y = int(obj_x), int(obj_y)
# data = data[y - hsize : y + hsize, x - hsize : x + hsize]

median = np.median(data)
std = np.median(np.abs(data - median))

fig = plt.figure()
ax = fig.add_subplot(111)
tr = transforms.Affine2D().rotate_deg(220)
ax.imshow(
    data, vmax=median + 7 * std, vmin=median - 3 * std, origin="lower", cmap="gray"
)
plt.plot(obj_x, obj_y, "oy")

csv_file = os.path.join(base_path, "GAIA_catalog.csv")
ss = pd.read_csv(csv_file)
for ra, dec in zip(ss["ra"], ss["dec"]):
    x, y, mjd = get_obj_coords(fits_file, ra, dec, unit=u.deg)
    plt.plot(x, y, "r+")


plt.xlabel("X axis")
plt.ylabel("Y axis")
plt.show()
