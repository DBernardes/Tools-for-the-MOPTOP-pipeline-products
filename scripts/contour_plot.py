#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import get_coords_in_series
import numpy as np
import pandas as pd
from scipy import stats
from sklearn import linear_model

star_name = "BD+32 3739"
experiment = "several positions in image/20230910"
camera = 4
alpha = 0.7
fontsize = 9
base_path = os.path.join(
    "..", "..", "Pol charact MOPTOP", "Low polarized stars", star_name, experiment
)


def prepare_data(new_path, parameter, filter="V"):
    df = pd.read_csv(new_path)
    rows = df.loc[df["wave"] == f"MOP-{filter}"]
    rows = rows.sort_values(["x", "y"], axis=0)

    x = np.asanyarray(rows["x"]).reshape(4, 4)
    y = np.asanyarray(rows["y"]).reshape(4, 4)
    val = np.asanyarray(rows[f"{parameter}_avg"]).reshape(4, 4)

    return x, y, val


def fit_plane(x, y, z):
    x1, y1, z1 = x.flatten(), y.flatten(), z.flatten()
    X_data = np.array([x1, y1]).reshape((-1, 2))
    Y_data = z1
    reg = linear_model.LinearRegression().fit(X_data, Y_data)
    a, b = reg.coef_
    c = reg.intercept_
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y + c
    return X, Y, Z


def calc_spearman(x, y, val):
    res = stats.spearmanr((x, y), val, axis=1)
    coor_x, coor_y, _ = res.statistic[-1]
    pval_x, pval_y, _ = res.pvalue[-1]

    return (
        res,
        pval_x,
        pval_y,
        coor_x,
        coor_y,
    )


def plot_data(ax, x, y, val):
    n = 5
    median = np.median(val)
    std = np.median(np.abs(val - median))
    levels = np.linspace(median - 3 * std, median + n * std, 15)

    cs = ax.contourf(x, y, val, levels=levels, alpha=0.75)
    # ax.contour(cs, colors="k", ls="-")
    ax.grid(c="k", ls="-", alpha=0.3)
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_title(f"{['q', 'u'][idx]} values")
    fig.colorbar(cs)

    return


import matplotlib.pyplot as plt
import numpy as np

# Data to plot.
new_path = os.path.join(base_path, "reduced", star_name, "manipulated_data.csv")

fig, axs = plt.subplots(1, 2)
for idx, parameter in enumerate(["q", "u"]):
    ax = axs[idx]
    x, y, val = prepare_data(new_path, parameter, "I")
    plot_data(ax, x, y, val)


plt.show()
