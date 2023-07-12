from tools import track_obj_over_images, _get_coords_in_series
import os
import matplotlib.pyplot as plt

min, max =  59685, 59705
star_name = 'BD+32 3739'
camera = 4
base_path = os.path.join('..', '..', 'Low polarized stars', star_name)
new_path  = os.path.join(base_path, star_name, 'selected_files', f'{min}-{max}', f'cam{camera}', f'I Filter')
# xcoord, ycoord, mjd = track_obj_over_images(new_path)

# plt.plot(mjd, xcoord, 'o-', label='x')
# plt.plot(mjd, ycoord, 'o-', label='y')
# plt.show()

date = _get_coords_in_series(new_path, ['2022-04-22'], 59691.148004)
print(date)