#!/usr/bin/env python

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




from tools import sort_qu_per_filter, low_polarized_stars, sigma_clipping
import os
import numpy as np
import pandas as pd


mean_qu_values = {'star':[], 'filter':[],'q':[], 'std_q':[], 'u':[], 'std_u':[]}

for star in low_polarized_stars.keys():
    base_path = os.path.join('..', '..', 'Pol charact MOPTOP', 'Low polarized stars', star, 'reduced', star )
    csv_file_name = os.path.join(base_path, 'manipulated_data.csv')
    qu_dict = sort_qu_per_filter(csv_file_name)
    for filter in qu_dict.keys():
        q = qu_dict[filter]['q']
        u = qu_dict[filter]['u']
        q, u = sigma_clipping(q, u)

        meanq, meanu = np.mean(q), np.mean(u)
        stdq, stdu = np.std(q), np.std(u)
        mean_qu_values['star'].append(star)
        mean_qu_values['filter'].append(filter)
        mean_qu_values['q'].append(meanq)
        mean_qu_values['std_q'].append(stdq)
        mean_qu_values['u'].append(meanu)
        mean_qu_values['std_u'].append(stdu)  

csv_file = os.path.join(base_path, '..', '..', '..', 'mean_qu_values.csv')
df = pd.DataFrame.from_dict(mean_qu_values)
df = df[df['q'].notna()]
pd.DataFrame.to_csv(df, csv_file, index=False)

