import os
import pandas as pd
from photometry import Photometry
from scripts.tools import sort_files
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from photometry import Photometry
from fits_files import FITS_files_manager
from math import log10
import sbpy
import collections

star_name = "GRB 230818A"
experiment = "first set"
object_type = "Scientific objects"
base_path = os.path.join(
    "..", "Pol charact MOPTOP", object_type, star_name, experiment, "reduced", star_name
)
# camera = 3
# csv_file = os.path.join(base_path, "..", "setup", f"objects coordinates copy.csv")
# objects = pd.read_csv(csv_file)
# objects = objects.loc[objects["name"] == "candidate1"]
# objects = objects.loc[:, ["name", f"ra_cam{camera}", f"dec_cam{camera}"]]
# objects = pd.DataFrame.from_dict(objects)
# file = f"4_e_20230818_5_10_16_1.fits"
# file_path = os.path.join(base_path, file)
# phot = Photometry(file_path, objects)
# phot.reset_object_coords()
# # phot.calculate_star_radius(coeff_radius_fwhm=2)
# phot.star_radius = 4.16
# phot.calc_sky_photons()
# phot.calc_star_photons()

# for obj in phot.obj_list:
#     print(obj)
#     print(phot.star_radius)

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
