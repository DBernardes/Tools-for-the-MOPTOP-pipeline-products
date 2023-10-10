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

star_name = "GRB 1163401"
experiment = "first set"
object_type = "Scientific objects"
base_path = os.path.join(
    "..", "Pol charact MOPTOP", object_type, star_name, experiment, "combined"
)

camera = 3
csv_file = os.path.join(
    base_path, "..", "setup", f"objects coordinates_combined images.csv"
)
objects = pd.read_csv(csv_file)

objects = objects.loc[objects["name"] == "candidate11"]
file = "3_e_20230408_10_6_16_1.fits"
file_path = os.path.join(base_path, file)
phot = Photometry(file_path, objects)
phot.reset_object_coords()
phot.calc_psf_radius()
phot.calc_sky_photons()
phot.calc_psf_photons()
for obj in phot.obj_list:
    print(obj)

# -----------------------------------------------------------------
