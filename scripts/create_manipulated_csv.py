#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os

from tools import manipulate_csv_file

star_name = "BD+32 3739"
obs_date = "20230909"
src_path = os.path.join("..", "..", "zpol stars", obs_date, star_name)
csv_file = os.path.join(src_path, "reduced_data.csv")
manipulate_csv_file(csv_file, os.path.join(src_path, ".."))
