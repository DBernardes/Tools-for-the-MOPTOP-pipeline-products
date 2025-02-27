#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn import linear_model
from tools import calculate_qu_low_standards, get_coords_in_series

star_name = "BD+32 3739"
base_path = os.path.join("..", "..", "zpol stars", "characterizations")


def test_plane_coefs(ax, parameter, filter):
    csv_path = os.path.join("csv", "plane_coefficients.csv")
    df = pd.read_csv(csv_path)
    df = df.loc[df["filter"] == filter]
    X = np.linspace(0, 1024, 10)
    Y = np.linspace(0, 1024, 10)
    Z = (
        df[parameter + "a"].values * X
        + df[parameter + "b"].values * Y
        + df[parameter + "c"].values
    )
    ax.plot(
        X,
        Y,
        Z,
        color="k",
        marker="o",
        alpha=0.5,
    )
    Y = np.linspace(1024, 0, 10)
    Z = (
        df[parameter + "a"].values * X
        + df[parameter + "b"].values * Y
        + df[parameter + "c"].values
    )
    ax.plot(
        X,
        Y,
        Z,
        color="k",
        marker="o",
        alpha=0.5,
    )
    return


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
    rows = rows[1:]
    (x, y, val, val_err) = (
        np.asanyarray(rows["x_pix"]),
        np.asanyarray(rows["y_pix"]),
        np.asanyarray(rows[f"{parameter}_avg"]),
        np.asanyarray(rows[f"{parameter}_err"]),
    )

    return x, y, val - lit_param, val_err


def fit_plane(x, y, z):
    x1, y1, z1 = x.flatten(), y.flatten(), z.flatten()
    X_data = np.array([x1, y1]).transpose()
    Y_data = z1
    reg = linear_model.LinearRegression().fit(X_data, Y_data)
    a, b = reg.coef_
    c = reg.intercept_
    X, Y = np.meshgrid(x, y)
    Z = a * X + b * Y + c
    return X, Y, Z, a, b, c


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


def plot_data(ax, x, y, val, val_err, parameter, filter):
    color = "r"
    if parameter == "q":
        color = "b"
        ax.set_title(f"Filter {filter}")

    ax.errorbar(
        x,
        y,
        val,
        val_err,
        color=color,
        marker="o",
        alpha=0.5,
        # label=f"{parameter}, corr:({coor_x:.3f},{coor_y:.3f}), pval:({pval_x:.2e},{pval_y:.2e})",
    )
    ax.set_xlim(0, 1024)
    ax.set_ylim(0, 1024)
    ax.set_xlabel("X (pix)")
    ax.set_ylabel("Y (pix)")
    # ax.set_zlabel(f"\n{parameter} values")
    ax.invert_yaxis()


fig = plt.figure()
for idx2, parameter in enumerate(["q", "u"]):
    for idx, filter in enumerate(["B", "V", "R", "I", "L"]):
        if filter in ["B", "L"]:
            star_name = "GD319"
        else:
            star_name = "BD+32 3739"
        new_path = os.path.join(base_path, f"filter {filter}.csv")
        x, y, val, val_err = prepare_data(new_path, parameter, filter, star_name)

        ax = fig.add_subplot(2, 5, idx2 * 5 + idx + 1, projection="3d")
        plot_data(ax, x, y, val, val_err, parameter, filter)
        # test_plane_coefs(ax, parameter, filter)
        X, Y, Z, a, b, c = fit_plane(x, y, val)
        ax.plot_surface(
            X,
            Y,
            Z,
            color=["b", "r"][idx2],
            alpha=0.1,
        )
        if parameter == "q":
            plt.title(
                f"Filter {filter}\n\nA={a:9.2e}\nB={b:9.2e}\nC={c:9.2e}", fontsize=10
            )
        else:
            plt.title(f"A={a:9.2e}\nB={b:9.2e}\nC={c:9.2e}", fontsize=10)
        if idx == 0:
            ax.zaxis.set_rotate_label(False)
            ax.set_zlabel(f"{parameter} values", rotation=90)

file = os.path.join(base_path, "3D_plot.png")
# plt.savefig(file, dpi=300)
plt.show()
