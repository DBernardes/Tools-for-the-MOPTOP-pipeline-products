import os
import matplotlib.pyplot as plt
from tools import calculate_mean_images,calculate_maximum_images


star_name = 'BD+32 3739'
src_path = os.path.join('..', 'Low polarized stars', star_name, star_name, 'selected_files', 'cam4')


for filter in ['V', 'R', 'I']:
    new_path = os.path.join(src_path, f'{filter} Filter')
    mean = calculate_mean_images(new_path)
    plt.plot(mean, 'o-', alpha = 0.5, label=f'{filter}')

plt.xlabel('Image number')
plt.ylabel('Counts level (ADU)')
plt.title('Mean of the pixels')
plt.legend()
plt.show()