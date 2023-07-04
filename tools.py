#!/usr/bin/env python

"""tools.py: this file has a set of functions developed for the analysis of the products obtained running the MOPTOP pipeline in a data set."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




import os
import pandas as pd
from math import sqrt, atan
import numpy as np
from copy import copy



def manipulate_csv_file(path):
    """Manipulate the raw_data csv file

    This function saves another csv file.
    This file has onyl the relevant information found in the raw_data file

    Args:
        path (str): path of the raw_data file
    """
    df = pd.read_csv(path)
    df = df[df['q_avg'].notna()]
    for column_name in df.columns.values:
        if column_name not in ['q_avg', 'q_err','u_avg','u_err', 'date', 'mjd', 'wave']:
            df = df.drop(column_name, axis=1)
    base_path = path.split('/')[:-1]
    new_path = os.path.join(*base_path, 'manipulated_data.csv')
    pd.DataFrame.to_csv(df, new_path, index=False)


def calculate_polarization(path):
    """Calculate the polarization of an object

    Given the q and u Stokes parameter in a csv file,
    this function calculates the percentage and the orientation of the polarization
    of an astronomical object.

    Args:
        path (str): path of the csv file
    """
    pol_dict = {}
    for filter in ['B', 'V', 'R', 'I', 'L']:
        pol_dict[filter] = {'mjd':[],
            'pol':[],
            'std':[],
            'pol_angle':[]}

    df = pd.read_csv(path)
    for _, row in df.iterrows():
        filter = row['wave'][-1]
        mjd = row['mjd']
        q = row['q_avg']
        q_err = row['q_err']
        u = row['u_avg']
        u_err = row['u_err']
        pol = sqrt(q**2 + u**2)
        std = sqrt((q/pol)**2 * q_err**2 + (u/pol)**2 * u_err)
        pol_angle = np.rad2deg(0.5 * atan(u/q))
        pol_dict[filter]['mjd'].append(mjd)
        pol_dict[filter]['pol'].append(pol)
        pol_dict[filter]['std'].append(std)
        pol_dict[filter]['pol_angle'].append(pol_angle)
    return pol_dict