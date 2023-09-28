#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


from tools import (
    get_instrumental_polarization,
    calculate_polarization_and_phase,
    novel_pol_error,
)
import os
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import pandas as pd


star_name = "GRB 230818A"
experiment = "all data/first set/combined images/2 positions"
object_type = "Scientific objects"
base_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    object_type,
    star_name,
    experiment,
    "reduced",
    star_name,
)
filter = "R"
fig, axs = plt.subplots(2, 1, figsize=(18, 5), sharex=True)
axs[0].set_title(f"Filter {filter}")


csv_file_name = os.path.join(base_path, "manipulated_data.csv")
df = pd.read_csv(csv_file_name)
q, u = df["q_avg"], df["u_avg"]
q_err, u_err = df["q_err"], df["u_err"]
x, y = df["x"], df["y"]
mjd = df["mjd"]
mjd = (mjd - mjd[0]) * 24 * 60


q_inst = get_instrumental_polarization(x, y, filter, "q")
u_inst = get_instrumental_polarization(x, y, filter, "u")
q_inst_err = np.std(q_inst)
u_inst_err = np.std(u_inst)
q -= q_inst
u -= u_inst
q_err = np.sqrt(q_err**2 + q_inst_err**2)
u_err = np.sqrt(u_err**2 + u_inst_err**2)
pol, pol_err, phase, phase_err = calculate_polarization_and_phase(q, q_err, u, u_err)
pol, pol_err, p_mas, p_mas_err, p_min, p_max = novel_pol_error(q, q_err, u, u_err)


ax = axs[0]
ax.errorbar(mjd, pol, pol_err, fmt="bo", alpha=0.5, label="naive estimator")
ax.errorbar(mjd, p_mas, p_mas_err, fmt="ro", alpha=0.5, label="mas estimator")
ax.set_ylabel(f"Polarization (%)")
ax.legend()
ax = axs[1]
ax.plot(mjd, phase, "bo", alpha=0.5)
ax.set_ylabel("Phase (deg)")
ax.set_xlabel("Time (min)")


plt.show()
