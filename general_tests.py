import os
import pandas as pd
from photometry import Photometry
from scripts.tools import sort_files
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from photometry import Photometry

star_name = "GRB 230818A"
experiment = "all data"
_set = "first"
camera = 3
src_path = os.path.join(
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
    experiment,
    star_name,
    f"{_set} set",
)


csv_file = os.path.join(src_path, "..", f"objects coordinates.csv")
df = pd.read_csv(csv_file)
objects = {
    "name": df["name"],
    "ra": df[f"ra_{_set}_set_cam{camera}"],
    "dec": df[f"dec_{_set}_set_cam{camera}"],
}
objects = pd.DataFrame.from_dict(objects)

objects_photometry = {}
for obj_name in objects["name"]:
    objects_photometry[obj_name] = {
        "mjd": [],
        "xcoord": [],
        "ycoord": [],
        "star_photons": [],
    }

file = "3_e_20230818_5_3_1_1.fits"
file_path = os.path.join(src_path, file)

phot = Photometry(file_path, objects, 20)
phot.reset_object_coords()
phot.calc_psf_radius()
phot.calc_sky_photons()
phot.calc_psf_photons()
mjd = phot.get_mjd()

for _object in phot.obj_list:
    print(repr(_object))
