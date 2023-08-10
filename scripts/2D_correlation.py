#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
from tools import get_coords_in_series
import numpy as np
import pandas as pd
from scipy import stats

star_name = "BD+32 3739"
experiment = "several positions in image/dense field"
camera = 4
alpha = 0.7
fontsize = 9
base_path = os.path.join(
    "..", "..", "Pol charact MOPTOP", "Low polarized stars", star_name, experiment
)


def calc_plot_parameters(df, x_str):
    x = df[x_str]
    q = df["q_avg"]
    u = df["u_avg"]
    q_err = df["q_err"]
    u_err = df["u_err"]
    res_q = stats.spearmanr(x, q)
    res_u = stats.spearmanr(x, u)
    return x, q, q_err, u, u_err, res_q, res_u


for idx, filter in enumerate(["V"]):
    new_path = os.path.join(base_path, "reduced", star_name, "manipulated_data.csv")
    df = pd.read_csv(new_path)
    rows = df.loc[df["wave"] == f"MOP-{filter}"]
    tmp = np.mean(rows["q_avg"] - rows["u_avg"])
    (
        x,
        y,
        q,
        u,
    ) = (
        np.asanyarray(rows["x"]),
        np.asanyarray(rows["y"]),
        np.asanyarray(rows["q_avg"]),
        np.asanyarray(rows["u_avg"]),
    )
    res_q = stats.spearmanr((x, y), q, axis=1)
    res_u = stats.spearmanr((x, y), u, axis=1)
    coor_qx, coor_qy, _ = res_q.statistic[-1]
    coor_ux, coor_uy, _ = res_u.statistic[-1]
    pval_qx, pval_qy, _ = res_q.pvalue[-1]
    pval_ux, pval_uy, _ = res_u.pvalue[-1]

    ax = plt.figure().add_subplot(projection="3d")
    ax.plot(
        x,
        y,
        q,
        zdir="z",
        c="b",
        marker="o",
        alpha=0.5,
        label=f"q, corr:({coor_qx:.3f},{coor_qy:.3f}), pval:({pval_qx:.2e},{pval_qy:.2e})",
    )
    ax.plot(
        x,
        y,
        u + tmp,
        zdir="z",
        c="r",
        marker="o",
        alpha=0.5,
        label=f"u+{tmp:.2f}, corr:({coor_ux:.3f},{coor_uy:.3f}), pval:({pval_ux:.2e},{pval_uy:.2e})",
    )
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.legend()


plt.show()
