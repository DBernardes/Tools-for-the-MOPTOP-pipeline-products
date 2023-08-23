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
src_path = os.path.join(
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
    experiment,
    star_name,
)

GOOD_IMAGE = "3_e_20230818_5_16_2_1.fits"
file = os.path.join(src_path, GOOD_IMAGE)
image = fits.getdata(file)
median = np.median(image)
std = np.median(np.abs(image - median))
plt.imshow(
    image, vmax=median + 7 * std, vmin=median - 3 * std, origin="lower", cmap="gray"
)

objects = [
    ("original", "19:03:33.2", "40:53:16.5"),
    ("comparison_1", "19:03:44.77", "40:52:07.4"),
    ("comparison_2", "19:03:42.6", "40:49:39"),
    ("candidate_1", "19:03:32", "40:53:11"),
    ("candidate_2", "19:03:33.39", "40:52:56.37"),
    ("candidate_3", "19:03:31", "40:53:34"),
]

phot = Photometry(file, objects)
for idx, _object in enumerate(phot.obj_list):
    name, x, y, *_ = _object.get_info()
    color = "b"
    if name == "original":
        color = "r"
    plt.plot(x, y, f"{color}o", alpha=0.25)
    plt.annotate(f"{idx+1}", (x * 0.99, y * 1.05), ha="right", va="bottom", fontsize=9)

plt.show()
