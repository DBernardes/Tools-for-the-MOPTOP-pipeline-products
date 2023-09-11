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
experiment = "all data"
_set = "first"
base_path = os.path.join(
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
    experiment,
    f"{_set} set",
    "combined images",
)

# ffiles = FITS_files_manager(base_path)
# shifts_file = os.path.join(base_path, "..", "shifts", "star_coords.csv")
# dest_path = os.path.join(base_path, "..", "combined images")

# ffiles.combine_images_by_run(dest_path, shifts_file)
# -----------------------------------------------------------------

camera = 3
csv_file = os.path.join(base_path, "..", f"objects coordinates.csv")
df = pd.read_csv(csv_file)
objects = {
    "name": ["object"],
    "ra": ["19:03:31.3029"],
    "dec": ["+40:53:43.022"],
}
objects = pd.DataFrame.from_dict(objects)

file = "3_e_run1.fits"
file_path = os.path.join(base_path, file)
phot = Photometry(file_path, objects, 20)
phot.reset_object_coords()
phot.calc_psf_radius()
phot.calc_sky_photons()
phot.calc_psf_photons()
for obj in phot.obj_list:
    print(obj)
