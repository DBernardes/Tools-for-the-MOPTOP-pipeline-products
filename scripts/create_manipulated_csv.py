#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


from tools import manipulate_csv_file
import os

star_name = "GRB 1163401"
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


csv_file_name = os.path.join(base_path, "raw_data.csv")
manipulate_csv_file(csv_file_name)
