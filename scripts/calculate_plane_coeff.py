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
    return a, b, c


for idx, filter in enumerate(["B", "V", "R", "I", "L"]):
    if filter in ["B", "L"]:
        star_name = "GD319"
    else:
        star_name = "BD+32 3739"
    print(filter)
    for idx2, parameter in enumerate(["q", "u"]):
        new_path = os.path.join(base_path, f"filter {filter}.csv")
        x, y, val, val_err = prepare_data(new_path, parameter, filter, star_name)
        a, b, c = fit_plane(x, y, val)
        print(a, b, c, sep=",")
