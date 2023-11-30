#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import numpy as np
import pandas as pd
import matplotlib.tri as tri
import matplotlib.cm as cm
from matplotlib import colormaps
from tools import calculate_qu_low_standards
import matplotlib.pyplot as plt
import numpy as np

alpha = 0.7
fontsize = 9
star_name = "HD14069"
experiment = "several positions in image/20231117/"
base_path = os.path.join(
    "..", "..", "Pol charact MOPTOP", "Low polarized stars", "csv contour plot"
)


def prepare_data(new_path, parameter, filter, star_name):
    lit_param = calculate_qu_low_standards(star_name, filter)
    if parameter == "q":
        lit_param = lit_param[0]
    elif parameter == "u":
        lit_param = lit_param[1]
    else:
        raise ValueError

    df = pd.read_csv(new_path)
    rows = df.loc[df["wave"] == f"MOP-{filter}"]
    rows = rows.sort_values(["x"], axis=0)

    x = np.asanyarray(rows["x"])
    y = np.asanyarray(rows["y"])
    val = np.asanyarray(rows[f"{parameter}_avg"])

    return x, y, val - lit_param


def plot_data(ax, x, y, val):
    ax.grid(c="k", ls="-", alpha=0.5)
    ax.set_xlim(0, 1024)
    ax.set_ylim(0, 1024)
    median = np.median(val)
    std = np.median(np.abs(val - median))
    cmap = colormaps["Blues"]
    levels = np.linspace(median - 3 * std, median + 5 * std, 15)
    triang = tri.Triangulation(x, y)
    # ax.tricontour(triang, val, colors="k", levels=levels, linestyles="solid", alpha=0.5)
    tcf = ax.tricontourf(triang, val, levels=levels, cmap=cmap)
    fig.colorbar(tcf)

    return


fig, axs = plt.subplots(1, 2, figsize=(12, 8), sharey="row", sharex="col")
axs[0].set_title("q values")
axs[1].set_title("u values")
axs[0].set_xlabel("X axis (pixels)")
axs[1].set_xlabel("X axis (pixels)")

for idx, filter in enumerate(["L"]):
    for idx2, parameter in enumerate(["q", "u"]):
        ax = axs[idx2]
        new_path = os.path.join(base_path, f"filter {filter}.csv")
        x, y, val = prepare_data(new_path, parameter, filter, star_name)
        plot_data(ax, x, y, val)
        if parameter == "q":
            ax.set_ylabel(f"Y axis (pixels)\n{filter} Filter")

file = os.path.join(base_path, "contour_plot.png")
plt.savefig(file, dpi=300)
plt.show()
