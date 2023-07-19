
from tools import sort_qu_per_filter, low_polarized_stars, high_polarized_stars, sigma_clipping
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sys import exit


star = 'HD14069'
base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star, 'reduced', star )
csv_file_name = os.path.join(base_path, 'manipulated_data.csv')

qu_dict = sort_qu_per_filter(csv_file_name)

q = qu_dict['B']['q']
u = qu_dict['B']['u']
q, u = sigma_clipping(q, u, 5, 3)
meanq, meanu = np.mean(q), np.mean(u)
plt.plot(q, u, 'o', alpha=0.25, label=f'{filter}')
plt.axhline(0, color='r', linestyle='--', alpha=0.25)
plt.axvline(0, color='r', linestyle='--', alpha=0.25)
plt.plot(meanq, meanu, '*')
plt.annotate(f'({meanq:.3f},{meanu:.3f})', (meanq*0.95,meanu), fontsize=10, ha='right')
plt.grid()
plt.show()