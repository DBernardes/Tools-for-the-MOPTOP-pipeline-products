#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


from tools import sort_qu_per_filter
import os
import matplotlib.pyplot as plt
import numpy as np
from sys import exit


star_name = "GRB 230818A"
experiment = "all data/first set"
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
fig, ax = plt.subplots(1, 1, figsize=(18, 5), sharex=True)
ax.set_title(f"Without combining the images")

for idx, filter in enumerate(["R"]):
    # ax = axs
    # ax.set_title(f'Filter {filter}')
    csv_file_name = os.path.join(base_path, "manipulated_data.csv")
    qu_dict = sort_qu_per_filter(csv_file_name)
    q = np.asarray(qu_dict[filter]["q"])
    u = np.asarray(qu_dict[filter]["u"])
    q_err = np.asarray(qu_dict[filter]["q_err"])
    u_err = np.asarray(qu_dict[filter]["u_err"])
    ax.errorbar(qu_dict[filter]["mjd"], q, yerr=q_err, fmt="bo", alpha=0.5, label=f"q")
    ax.errorbar(qu_dict[filter]["mjd"], u, yerr=u_err, fmt="ro", alpha=0.5, label=f"u")
    # ax.set_ylim(-0.1, 0.1)
    # ax.set_xlim(59600, 60100)
    ax.legend()
    ax.set_ylabel(f"Filter {filter}")
plt.xlabel("Time (MJD)")
plt.show()
