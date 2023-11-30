#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


from tools import manipulate_csv_file
import os

star_name = "HD14069"
experiment = "several positions in image/20231117/"
src_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    "Low polarized stars",
    star_name,
    experiment,
    "reduced",
    star_name,
)
csv_file = os.path.join(src_path, "raw_data.csv")
dest_path = os.path.join(src_path, "..", "..", "polarimetry")
manipulate_csv_file(csv_file, src_path)
