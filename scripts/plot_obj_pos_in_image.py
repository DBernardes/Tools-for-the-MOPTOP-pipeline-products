#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np
from tools import track_obj_over_images
import matplotlib.patches as patches

star_name = "HD14069"
experiment = "several positions in image/20231107"
src_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    "Low polarized stars",
    star_name,
    experiment,
    star_name,
)
fig, axs = plt.subplots(1, 2, figsize=(15, 7), sharex=True, sharey=True)
axs[0].set_ylabel("Y axis - camera 1")
axs[1].set_ylabel("Y axis - camera 2")

for camera in [3, 4]:
    for idx, filter in enumerate(["V"]):
        ax = axs[camera - 3]
        new_path = src_path
        xcoord, ycoord, _ = track_obj_over_images(new_path, f"{camera}_e")
        ax.plot(xcoord, ycoord, "o-", alpha=0.8)

        medx, medy = np.median(xcoord), np.median(ycoord)

        ax.axvline(512, color="r", linestyle="--", alpha=0.75)
        ax.axhline(512, color="r", linestyle="--", alpha=0.75)
        # Create a Rectangle patch
        rect = patches.Rectangle(
            (0, 0),
            1024,
            1024,
            linewidth=1,
            edgecolor="k",
            facecolor="none",
        )
        ax.add_patch(rect)
        if camera == 3:
            ax.set_title(f"Filter {filter}")
        else:
            ax.set_xlabel("X axis")


plt.xlim(0, 1024)
plt.ylim(0, 1024)
plt.show()
