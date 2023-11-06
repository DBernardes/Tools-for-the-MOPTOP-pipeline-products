#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import pandas as pd
from tools import calculate_polarization_and_phase


star_name = "GRB 230818A"
experiment = "first set"
object_type = "Scientific objects"
base_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    object_type,
    star_name,
    experiment,
    "polarization",
)
filter = "R"


csv_file_name = os.path.join(base_path, "cand2.csv")
df = pd.read_csv(csv_file_name)
df = df[:-1]
q, u = df["q_avg"], df["u_avg"]
q_err, u_err = df["q_err"], df["u_err"]
# pol, pol_err, *_ = calculate_polarization_and_phase(q, q_err, u, u_err)
x, y = df["x"], df["y"]
mjd = df["mjd"]
mjd = (mjd - mjd[0]) * 24 * 60

plt.errorbar(mjd, q, q_err, fmt="bo", alpha=0.5, label="q")
plt.errorbar(mjd, u, u_err, fmt="ro", alpha=0.5, label="u")
# plt.errorbar(mjd, pol, pol_err, fmt="go", alpha=0.5, label="pol")
plt.ylabel(f"Stokes parameters")
plt.xlabel("Time (min)")
plt.legend()
plt.show()
