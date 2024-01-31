#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import pandas as pd
from matplotlib import colormaps
from tools import calculate_qu_low_standards

alpha = 0.7
fontsize = 9
star_name = "GD319"
base_path = os.path.join("..", "..", "zpol stars", "characterizations")


def prepare_data(new_path, parameter, filter, star_name):
    lit_param = calculate_qu_low_standards(star_name, filter)
    if parameter == "q":
        lit_param = lit_param[0]
    elif parameter == "u":
        lit_param = lit_param[1]
    else:
        raise ValueError

    df = pd.read_csv(new_path)
    rows = df.loc[df["wave"] == f"{filter}"]
    # rows = rows.sort_values(["x_pix"], axis=0)

    x = np.asanyarray(rows["x_pix"])
    y = np.asanyarray(rows["y_pix"])
    val = np.asanyarray(rows[f"{parameter}_avg"])
    # x, y, val = zip(*(x, y, val))

    return x, y, val - lit_param


def plot_data(ax, x, y, val):
    ax.grid(c="k", ls="-", alpha=0.5)
    ax.set_xlim(0, 1024)
    ax.set_ylim(0, 1024)
    median = np.median(val)
    std = np.std(val)
    cmap = colormaps["Blues"]
    levels = np.linspace(median - 4 * std, median + 4 * std, 15)
    triang = tri.Triangulation(x, y)
    # ax.tricontour(triang, val, colors="k", levels=levels, linestyles="solid", alpha=0.5)
    tcf = ax.tricontourf(triang, val, levels=levels, cmap=cmap)
    fig.colorbar(tcf)

    return


fig, axs = plt.subplots(4, 2, figsize=(12, 8), sharey="row", sharex="col")
axs[0, 0].set_title("q values")
axs[0, 1].set_title("u values")
axs[3, 1].set_xlabel("X axis (pixels)")
axs[3, 0].set_xlabel("X axis (pixels)")

for idx, filter in enumerate(["B", "V", "R", "L"]):
    for idx2, parameter in enumerate(["q", "u"]):
        ax = axs[idx, idx2]
        new_path = os.path.join(base_path, f"filter {filter}.csv")
        x, y, val = prepare_data(new_path, parameter, filter, star_name)
        plot_data(ax, x, y, val)
        if parameter == "q":
            ax.set_ylabel(f"Y axis (pixels)\n{filter} Filter")

file = os.path.join(base_path, "contour_plot.png")
plt.savefig(file, dpi=300)
plt.show()
