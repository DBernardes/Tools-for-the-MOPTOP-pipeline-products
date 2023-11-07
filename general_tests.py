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

star_name = "GRB 230818A"
experiment = "first set"
object_type = "Scientific objects"
base_path = os.path.join(
    "..", "Pol charact MOPTOP", object_type, star_name, experiment, star_name
)

camera = 3
csv_file = os.path.join(base_path, "..", "setup", f"objects coordinates.csv")
objects = pd.read_csv(csv_file)
objects = objects.loc[objects["name"] == "candidate1"]
objects = objects.drop(["ra_cam4", "dec_cam4"], axis=1)
objects = pd.DataFrame.from_dict(objects)
file = "3_e_20230818_5_4_1_1.fits"
file_path = os.path.join(base_path, file)
phot = Photometry(file_path, objects)
phot.reset_object_coords()
phot.calculate_star_radius(coeff_radius_fwhm=2)

for obj in phot.obj_list:
    print(obj)
    print(phot.star_radius)

# -----------------------------------------------------------------


# ffile_man = FITS_files_manager(base_path)
# dest_path = os.path.join(base_path, "..", "tmp")
# ffile_man.combine_images_by_rotor_position(dest_path, 2)
