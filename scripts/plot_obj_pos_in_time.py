#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import track_obj_over_images
import numpy as np
import pandas as pd

min, max = 59775, 60051
star_name = "GRB 230818A"
experiment = "all data"
base_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
    experiment,
    star_name,
)
fig, axs = plt.subplots(1, 2, figsize=(15, 7))

for camera in [3, 4]:
    new_path = base_path
    xcoord, ycoord, mjd = track_obj_over_images(new_path, f"{camera}_e")
    ax = axs[camera - 3]
    ax.plot(mjd, xcoord, "o", label=f"x coord")
    ax.plot(mjd, ycoord, "o", label=f"y coord")
    ax.set_ylabel("Pixels")
    ax.set_xlim(mjd[0], mjd[-1])
    ax.legend()

    # new_path = os.path.join(base_path, "reduced", star_name, "manipulated_data.csv")
    # df = pd.read_csv(new_path)
    # rows = df.loc[df["wave"] == f"MOP-{filter}"]
    # photons = rows[f"s1_src_cam{camera-2}"]
    # mean = np.median(photons)
    # std = np.std(np.abs(photons - mean))

    # ax = axs[1, idx]
    # ax.plot(rows["mjd"], photons, "o", label=f"cam{camera}")
    # ax.set_ylabel("Photons")
    # ax.set_xlim(mjd[0], mjd[-1])
    # ax.set_ylim(mean - 1.5 * std, mean + 1.5 * std)
    # ax.invert_yaxis()
    # ax.legend()

    # ax = axs[2, idx]
    # q = rows["q_avg"]
    # u = rows["u_avg"]
    # tmp = np.mean(q - u)
    # ax.plot(rows["mjd"], q, "o", label=f"q")
    # ax.plot(rows["mjd"], u + tmp, "o", label=f"u + {tmp:.2f}")
    # ax.set_xlabel("Time (MJD)")
    # ax.set_xlim(mjd[0], mjd[-1])
    # ax.legend()


plt.show()
