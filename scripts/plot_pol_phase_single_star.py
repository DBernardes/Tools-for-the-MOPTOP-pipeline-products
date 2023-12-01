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
experiment = "first set"
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
fig, axs = plt.subplots(2, 1, figsize=(6, 5), sharex=True)

csv_files = [file for file in os.listdir(base_path) if "cand" in file]

for file in csv_files:
    csv_file_name = os.path.join(base_path, file)
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
    pol, pol_err, phase, phase_err = calculate_polarization_and_phase(
        q, q_err, u, u_err
    )
    pol, pol_err, p_mas, p_mas_err, p_min, p_max = novel_pol_error(q, q_err, u, u_err)

    ax = axs[0]
    ax.errorbar(mjd, p_mas, p_mas_err, fmt="o", alpha=0.5, label=file.split(".")[0])
    ax.set_ylabel(f"Polarization (%)")
    ax.set_ylim(0)
    ax.legend()
    ax = axs[1]
    ax.errorbar(mjd, phase, phase_err, fmt="o", alpha=0.5, label=file.split(".")[0])
    ax.set_ylabel("Phase (deg)")
    ax.set_xlabel("Time (min)")
    ax.legend()

# -------------------------------------------------------------------------------------------

# dest_files_path = os.path.join(base_path, "polarimetry")
# csv_file = os.path.join(dest_files_path, "polarization.csv")

# df = pd.read_csv(csv_file)
# q, u = df["q"], df["u"]
# q_err, u_err = df["q_err"], df["u_err"]
# mjd = df["mjd"]
# x, y = df["x"], df["y"]
# mjd = (mjd - mjd[0]) * 24 * 60

# q_inst = get_instrumental_polarization(x, y, "R", "q")
# u_inst = get_instrumental_polarization(x, y, "R", "u")
# q_inst_err = np.std(q_inst)
# u_inst_err = np.std(u_inst)
# q -= q_inst
# u -= u_inst
# q_err = np.sqrt(q_err**2 + q_inst_err**2)
# u_err = np.sqrt(u_err**2 + u_inst_err**2)
# pol, pol_err, phase, phase_err = calculate_polarization_and_phase(q, q_err, u, u_err)
# pol, pol_err, p_mas, p_mas_err, p_min, p_max = novel_pol_error(q, q_err, u, u_err)
# ax = axs[0]
# ax.errorbar(mjd, p_mas, p_mas_err, fmt="ro", alpha=0.5, label="SPARC4")
# ax.legend()
# ax = axs[1]
# ax.errorbar(mjd, phase, phase_err, fmt="ro", alpha=0.5, label="SPARC4")
# ax.legend()


fig_file = os.path.join(base_path, "..", "..", "plots", "moptop_pol.png")
plt.savefig(fig_file, dpi=300)
plt.show()
