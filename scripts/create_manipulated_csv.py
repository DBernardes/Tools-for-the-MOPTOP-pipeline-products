#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


from tools import manipulate_csv_file
import os

star_name = "HD14069"
experiment = "several positions in image/20231107"
object_type = "Low polarized stars"
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

dest_path = os.path.join(base_path, *2 * [".."], "polarization")
csv_file_name = os.path.join(base_path, "raw_data.csv")
manipulate_csv_file(csv_file_name, base_path)
