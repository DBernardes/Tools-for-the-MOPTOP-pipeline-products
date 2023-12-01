#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


from tools import manipulate_csv_file
import os

star_name = "GRB 230818A"
experiment = "first set"
src_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
    experiment,
    "reduced",
    star_name,
)
csv_file = os.path.join(src_path, "raw_data.csv")
manipulate_csv_file(csv_file, src_path, "cand8")
