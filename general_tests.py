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
    "..", "Pol charact MOPTOP", object_type, star_name, experiment, "all data"
)


# ffiles = FITS_files_manager(base_path)
# shifts_file = os.path.join(base_path, "..", "setup", "star_coords.csv")
# dest_path = os.path.join(base_path, "..", "combined images", "2 positions", star_name)

# ffiles.combine_images_by_rotor_position(
#     dest_path, shifts_file, nruns=2, use_moptp_name=True
# )

# -----------------------------------------------------------------

camera = 3
csv_file = os.path.join(base_path, "..", "setup", f"objects coordinates.csv")
df = pd.read_csv(csv_file)
objects = {
    "name": ["object"],
    "ra": ["6:34:22.0457"],
    "dec": ["+49:48:36.379"],
}
objects = pd.DataFrame.from_dict(objects)
file = "3_e_20230116_19_1_1_1.fits"
file_path = os.path.join(base_path, file)
phot = Photometry(file_path, objects)
phot.reset_object_coords()
phot.calc_psf_radius()
for obj in phot.obj_list:
    print(obj)

# -----------------------------------------------------------------
