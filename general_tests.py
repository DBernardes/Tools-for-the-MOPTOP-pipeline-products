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

star_name = "GRB 1149293"
experiment = "first set"
object_type = "Scientific objects"
base_path = os.path.join(
    "..", "Pol charact MOPTOP", object_type, star_name, experiment, "combined"
)

camera = 3
# csv_file = os.path.join(base_path, "..", "setup", f"objects coordinates.csv")
# objects = pd.read_csv(csv_file)
# objects = objects.loc[objects["name"] == "candidate7"]
base_path = os.path.join("tests", "images")
ra, dec = "19:03:40.0436", "+40:50:28.191"
objects = {"name": [star_name], "ra_cam3": [ra], "dec_cam3": [dec]}
objects = pd.DataFrame.from_dict(objects)
file = "3_e_20230818_5_1_1_1.fits"
file_path = os.path.join(base_path, file)
phot = Photometry(file_path, objects)
phot.reset_object_coords()

for obj in phot.obj_list:
    print(obj)

# -----------------------------------------------------------------


# ffile_man = FITS_files_manager("./tests/images")
# ffile_man._verify_rotor_positions()
