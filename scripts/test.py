
from tools import sort_qu_per_filter, low_polarized_stars, high_polarized_stars, sigma_clipping
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sys import exit
from circle_fit import taubinSVD


star = 'VICyg12'
base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars')
csv_file_name = os.path.join(base_path, 'mean_qu_values.csv')

 
df = pd.read_csv(csv_file_name)
for star in df['star']:
    if star not in ['GD319', 'HD14069', 'BD+32 3739']:
        df.drop(df[df['star'] == star].index, inplace = True)

rows = df.loc[df['filter'] == 'B']

print(rows)