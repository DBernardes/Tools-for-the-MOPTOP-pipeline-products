from tools import _sort_files
import os

star_name = 'BD+32 3739'
src_path = os.path.join('..', '..', 'Low polarized stars', star_name, star_name, 'selected_files', 'cam3', 'R Filter')
files = _sort_files(src_path, '.fits')

for file in files:print(file)
 