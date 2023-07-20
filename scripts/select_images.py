#!/usr/bin/env python

"""run.py: this is a scrit to run the functions of the tools.py file."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import select_images_keyword_interval, delete_file_keyword_value, select_images_keyword_value
import os
from sys import exit

min, max =  59860,59932 #parecido com o intervalo da outra estrela
star_name = 'VICyg12'
src_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'High polarized stars', star_name, star_name, 'LTOps')
dest_path = os.path.join(src_path, '..', 'selected_files', f'{min}-{max}')
# select_images_keyword_interval(src_path, dest_path, 'MJD', min, max)


for filter in ['V', 'R', 'I']:
    src_path_1 = os.path.join(dest_path, 'cam4')
    dest_path_1 = os.path.join(src_path_1, f'{filter} Filter')
    select_images_keyword_value(src_path_1, dest_path_1, 'FILTER1', f'MOP-{filter}')

