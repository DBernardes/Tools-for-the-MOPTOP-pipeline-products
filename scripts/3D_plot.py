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
experiment = "several positions in image/20230830"
camera = 4
alpha = 0.7
fontsize = 9
base_path = os.path.join(
    "..", "..", "Pol charact MOPTOP", "Low polarized stars", star_name, experiment
)


def prepare_data(new_path, parameter, filter="V"):
    df = pd.read_csv(new_path)
    rows = df.loc[df["wave"] == f"MOP-{filter}"]
    rows = rows.drop(15)
    (x, y, val, val_err) = (
        np.asanyarray(rows["x"]),
        np.asanyarray(rows["y"]),
        np.asanyarray(rows[f"{parameter}_avg"]),
        np.asanyarray(rows[f"{parameter}_err"]),
    )

    return x, y, val, val_err


def fit_plane(x, y, z):
    x1, y1, z1 = x.flatten(), y.flatten(), z.flatten()
    X_data = np.array([x1, y1]).transpose()
    Y_data = z1
    reg = linear_model.LinearRegression().fit(X_data, Y_data)
    a, b = reg.coef_
    c = reg.intercept_
    X, Y = np.meshgrid(x, y)
    print(a, b, c)
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


def plot_data(ax, x, y, val, val_err, coor_x, coor_y, pval_x, pval_y, parameter):
    color = "b"
    if parameter == "u":
        color = "r"
    ax.errorbar(
        x,
        y,
        val,
        val_err,
        color=color,
        marker="o",
        alpha=0.5,
        label=f"{parameter}, corr:({coor_x:.3f},{coor_y:.3f}), pval:({pval_x:.2e},{pval_y:.2e})",
    )


fig = plt.figure(figsize=plt.figaspect(0.5))
new_path = os.path.join(base_path, "reduced", star_name, "manipulated_data.csv")
for idx, parameter in enumerate(["q", "u"]):
    x, y, val, val_err = prepare_data(new_path, parameter, "V")

    (
        res,
        pval_x,
        pval_y,
        coor_x,
        coor_y,
    ) = calc_spearman(x, y, val)
    ax = fig.add_subplot(1, 2, idx + 1, projection="3d")
    plot_data(ax, x, y, val, val_err, coor_x, coor_y, pval_x, pval_y, parameter)
    X, Y, Z = fit_plane(x, y, val)
    ax.plot_surface(
        X,
        Y,
        Z,
        color=["b", "r"][idx],
        alpha=0.1,
    )
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.legend()


plt.show()
