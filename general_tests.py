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
experiment = "all data/first set/"
object_type = "Scientific objects"
base_path = os.path.join(
    "..",
    "Pol charact MOPTOP",
    object_type,
    star_name,
    experiment,
    star_name,
)


ffiles = FITS_files_manager(base_path)
shifts_file = os.path.join(base_path, "..", "setup", "star_coords.csv")
dest_path = os.path.join(base_path, "..", "combined images", "2 positions", star_name)

ffiles.combine_images_by_rotor_position(
    dest_path, shifts_file, nruns=2, use_moptp_name=True
)

# -----------------------------------------------------------------

# camera = 3
# csv_file = os.path.join(base_path, "..", "setup", f"objects coordinates.csv")
# df = pd.read_csv(csv_file)
# objects = {
#     "name": ["object"],
#     "ra": ["19:03:33.1120"],
#     "dec": ["+40:54:08.919"],
# }
# objects = pd.DataFrame.from_dict(objects)

# file = "4_e_20230818_5_2_2_1.fits"
# file_path = os.path.join(base_path, file)
# phot = Photometry(file_path, objects, 20)
# phot.reset_object_coords()
# phot.calc_psf_radius()
# phot.calc_sky_photons()
# phot.calc_psf_photons()
# for obj in phot.obj_list:
#     print(obj)

# -----------------------------------------------------------------
