import collections
import os
from math import log10

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sbpy

from fits_files import FITS_files_manager
from photometry import Photometry
from scripts.tools import sort_files

star_name = "GRB 230818A"
base_path = os.path.join("..", "Dados GRB", star_name)

camera = 4
radius = []
for i in range(1, 17):
    csv_file = os.path.join(base_path, "..", "setup", f"objects coordinates.csv")
    objects = pd.read_csv(csv_file)
    objects = objects.loc[objects["name"] == "candidate1"]
    objects = objects.loc[:, ["name", f"ra_cam{camera}", f"dec_cam{camera}"]]
    objects = pd.DataFrame.from_dict(objects)
    file = f"{camera}_e_20230818_5_2_{i}_1.fits"
    file_path = os.path.join(base_path, file)
    phot = Photometry(file_path, objects)
    phot.reset_object_coords()
    phot.calculate_star_radius(coeff_radius_fwhm=2)
    phot.calc_sky_photons()
    phot.calc_star_photons()
    radius.append(phot.star_radius)

print(np.mean(radius))


# -----------------------------------------------------------------


# ffile_man = FITS_files_manager(base_path)
# dest_path = os.path.join(base_path, "..", "tmp")
# ffile_man.combine_images_by_rotor_position(dest_path, 2)


# ----------------------------------------------------------------
# camera = 4
# csv_file = os.path.join(base_path, "..", "setup", f"objects coordinates.csv")
# objects = pd.read_csv(csv_file)
# objects = objects.loc[objects["name"] == "candidate1"]
# objects = objects.drop(["ra_cam4", "dec_cam4"], axis=1)
# ffiles = FITS_files_manager(base_path)

# tmp = []
# run_files = ffiles.get_images_by_run(10)
# for file in run_files[f"cam{camera}"]:
#     file = os.path.join(base_path, file.name)
#     phot = Photometry(file, objects)
#     phot.reset_object_coords()
#     star_radius = phot.calculate_star_radius(coeff_radius_fwhm=2)
#     print(file)
#     print(star_radius)
