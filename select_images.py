#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import select_images_keyword_interval, delete_file_keyword_value, select_images_keyword_value
import os
from sys import exit

star_name = 'BD+32 3739'
base_path = os.path.join('..', 'Low polarized stars', star_name, star_name)
src_path = os.path.join(base_path, 'selected_files', 'cam4')
dest_path = os.path.join(src_path, 'V Filter')
select_images_keyword_value(src_path, dest_path, 'FILTER1', 'MOP-V')


