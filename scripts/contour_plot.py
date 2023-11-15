#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.tri as tri
import matplotlib.cm as cm
from matplotlib import colormaps

star_name = "BD+32 3739"
experiment = "several positions in image/20230910"
camera = 4
alpha = 0.7
fontsize = 9
base_path = os.path.join(
    "..", "..", "Pol charact MOPTOP", "Low polarized stars", "csv contour plot"
)


def prepare_data(new_path, parameter, filter):
    df = pd.read_csv(new_path)
    rows = df.loc[df["wave"] == f"MOP-{filter}"]
    rows = rows.sort_values(["x"], axis=0)

    x = np.asanyarray(rows["x"])
    y = np.asanyarray(rows["y"])
    val = np.asanyarray(rows[f"{parameter}_avg"])

    return x, y, val


def plot_data(ax, x, y, val):
    ax.grid(c="k", ls="-", alpha=0.5)
    ax.set_xlim(0, 1024)
    ax.set_ylim(0, 1024)
    median = np.median(val)
    std = np.median(np.abs(val - median))
    cmap = colormaps["Blues"]
    levels = np.linspace(median - 3 * std, median + 5 * std, 15)
    triang = tri.Triangulation(x, y)
    ax.tricontour(triang, val, colors="k", levels=levels, linestyles="solid", alpha=0.5)
    tcf = ax.tricontourf(triang, val, levels=levels, cmap=cmap)
    fig.colorbar(tcf)

    return


def plot_data_1(ax, x, y, val):
    n = 5
    median = np.median(val)
    std = np.median(np.abs(val - median))
    levels = np.linspace(median - 3 * std, median + n * std, 15)

    cs = ax.contourf(x, y, val, levels=levels, alpha=0.75)
    ax.grid(c="k", ls="-", alpha=0.3)

    fig.colorbar(cs)

    return


import matplotlib.pyplot as plt
import numpy as np

# Data to plot.

fig, axs = plt.subplots(3, 2, figsize=(12, 8), sharey="row", sharex="col")
axs[0, 0].set_title("q values")
axs[0, 1].set_title("u values")
axs[2, 0].set_xlabel("X axis (pixels)")
axs[2, 1].set_xlabel("X axis (pixels)")

for idx, filter in enumerate(["V", "R", "I"]):
    for idx2, parameter in enumerate(["q", "u"]):
        ax = axs[idx, idx2]
        new_path = os.path.join(base_path, f"filter {filter}.csv")
        x, y, val = prepare_data(new_path, parameter, filter)
        plot_data(ax, x, y, val)
        if parameter == "q":
            ax.set_ylabel(f"Y axis (pixels)\n{filter} Filter")

file = os.path.join(base_path, "contour_plot.png")
plt.savefig(file, dpi=300)
plt.show()
