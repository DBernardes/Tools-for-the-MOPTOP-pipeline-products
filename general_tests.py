import os
import pandas as pd
from photometry import Photometry
from scripts.tools import sort_files
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from photometry import Photometry
from fits_files import FITS_files_manager

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
    star_name,
)

ffiles = FITS_files_manager(base_path)
dest_path = os.path.join(base_path, "..", "combined images")
# for obj in ffiles.get_images_by_run(7)["cam4"]:
#     print(obj)

# -----------------------------------------------------------------

camera = 4
csv_file = os.path.join(base_path, "..", f"objects coordinates.csv")
df = pd.read_csv(csv_file)
objects = {
    "name": df["name"],
    "ra": df[f"ra_{_set}_set_cam{camera}"],
    "dec": df[f"dec_{_set}_set_cam{camera}"],
}
objects = pd.DataFrame.from_dict(objects)

for name in ffiles.get_images_by_run(7)["cam4"]:
    name = name.name
    file_path = os.path.join(base_path, name)
    phot = Photometry(file_path, objects, 20)
    phot.reset_object_coords()
    phot.calc_psf_radius()
    phot.calc_sky_photons()
    phot.calc_psf_photons()
    print(phot.obj_list[2])
