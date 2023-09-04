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
# dest_path = os.path.join(base_path, "..", "combined images")
# ffiles.combine_images_by_run(dest_path)


# -----------------------------------------------------------------

camera = 3
csv_file = os.path.join(base_path, "..", f"objects coordinates.csv")
df = pd.read_csv(csv_file)
objects = {
    "name": df["name"],
    "ra": df[f"ra_{_set}_set_cam{camera}"],
    "dec": df[f"dec_{_set}_set_cam{camera}"],
}
objects = pd.DataFrame.from_dict(objects)

file = "cam3_run14.fits"
file_path = os.path.join(base_path, file)
phot = Photometry(file_path, objects, 30)
phot.reset_object_coords()
phot.calc_psf_radius()
phot.calc_sky_photons()
phot.calc_psf_photons()
for obj in phot.obj_list:
    print(obj)
