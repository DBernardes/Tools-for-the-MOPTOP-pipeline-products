#!/usr/bin/env python

"""tools.py: this file has a set of functions developed for the analysis of the products obtained running the MOPTOP pipeline in a data set."""

__author__      = "Denis Bernardes"
__copyright__   = "Copyright 2023, Liverpool John Moores University"




import os
import pandas as pd
import numpy as np
import astropy.io.fits as fits
import shutil
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from stat import S_ISREG, ST_CTIME, ST_MODE
import glob
from numpy import ndarray
from sys import exit


low_palarized_stars = {'GD319':0.045, 'HD14069':0.111, 'BD+32 3739':0.039, 'HD212311':0.028, 'BD+28 4211':0.063, 'G191B2B':0.090, 'BD+33 2642':0.145}
#https://www.not.iac.es/instruments/turpol/std/zpstd.html

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
        if column_name not in ['q_avg', 'q_err','u_avg','u_err', 'date', 'mjd', 'wave', 's1_src_cam1', 's1_src_cam1_err', 's1_src_cam2', 's1_src_cam2_err', 'time']:
            df = df.drop(column_name, axis=1)
    df = df.sort_values(by=['wave', 'mjd'])
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
    qu_dict = sort_qu_per_filter(path)
    pol_dict = {'B':[], 'V':[], 'R':[], 'I':[], 'L':[]}
    for key in qu_dict.keys():
        mjd = np.asarray(qu_dict[key]['mjd'])
        q = np.asarray(qu_dict[key]['q'])
        u = np.asarray(qu_dict[key]['u'])
        pol = np.sqrt(np.square(q) + np.square(u))*100
        pol_dict[key].append(mjd)
        pol_dict[key].append(pol)
        
    return pol_dict

def sort_qu_per_filter(path:str)-> dict:
    """Sort the q and u values of the Stoke parameters for each filter

    Args:
        path (str): csv file path
    """
    df = pd.read_csv(path)
    qu_dict = {}
    for filter in ['B', 'V', 'R', 'I', 'L']:
        rows = df.loc[df['wave'] == f'MOP-{filter}']
        qu_dict[filter] = {'mjd': rows['mjd'],
            'q':rows['q_avg'],
            'u':rows['u_avg']
            }

    return qu_dict

def track_obj_over_images(path:str, tag:str='.fits'):
    """Track object over the images

    Args:
        path (str): directory of the images
    """
    files = _sort_files(path, tag)
    xcoord, ycoord, mjds = ndarray, ndarray, ndarray
    for file in files:
        x, y, mjd = get_obj_coords(path, file)
        xcoord.append(x)
        ycoord.append(y)
        mjds.append(mjd)
    return xcoord, ycoord, mjds

def get_obj_coords(path:str, file:str):
    hdr = fits.getheader(os.path.join(path, file))
    wcs = WCS(hdr)  
    coord = SkyCoord(f'{hdr["CAT-RA"]} {hdr["CAT-DEC"]}', frame='fk5', unit=(u.hourangle, u.deg))
    x, y = wcs.world_to_pixel(coord)
    return x, y, hdr['MJD']

def select_images_keyword_interval(path: str, dest_path: str, keyword:str, min:float = -np.infty, max:float = np.infty, tag:str='.fits')-> None:
    """Select those images in which the keyword is inside the min and max values range

    Args:
        src_path (str): path of the FITS files
        dest_path (str): destination of the selected files
        min (float, optional): minimum value for MJD. Defaults to None.
        max (float, optional): maximum value for MJD. Defaults to None.
    """
    files = _sort_files(path, tag)
    os.makedirs(dest_path, exist_ok=True)
    for file in files:
        src_file = os.path.join(path, file)
        hdr = fits.getheader(src_file)
        if min < hdr[keyword] < max:
            dest_file = os.path.join(dest_path, file)
            shutil.copyfile(src_file, dest_file)

def select_images_keyword_value(path: str, dest_path: str, keyword:str, value: float|int|str, tag:str='.fits')-> None:
    """Select those images in which the keyword matches the provided value

    Args:
        src_path (str): path of the FITS files
        dest_path (str): destination of the selected files
        keyword (str): header keyword
        value (float, int, str): value of the keyword
    """
    files = _sort_files(path, tag)
    os.makedirs(dest_path, exist_ok=True)
    for file in files:
        src_file = os.path.join(path, file)
        hdr = fits.getheader(src_file)
        if hdr[keyword] == value:
            dest_file = os.path.join(dest_path, file)
            shutil.copyfile(src_file, dest_file)

def delete_file_keyword_value(path: str, keyword:str, value: float | str | int, tag:str='.fits')-> None:
    """Delete files found in a folder based on the provided header keyword value

    Args:
        src_path (str): source path
        keyword (str): header keyword
        value (float, int, str): keyword value
    """
    files = _sort_files(path, tag)
    for file in files:
        src_file = os.path.join(path, file)
        hdr = fits.getheader(src_file)
        if hdr[keyword] == value:
            os.remove(src_file)
    return
 
def calculate_mean_images(path, tag:str='.fits'):
    files = _sort_files(path, tag)
    mean = []
    for file in files:
        image = fits.getdata(os.path.join(path, file))
        tmp = np.mean(image)
        mean.append(tmp)
    return np.asarray(mean)

def calculate_maximum_images(path, tag:str='.fits'):
    files = _sort_files(path, tag)
    _max = []
    for file in files:
        image = fits.getdata(os.path.join(path, file))
        tmp = np.max(image)
        _max.append(tmp)
    return np.asarray(_max)

def _sort_files(path:str, tag):
    files = [f for f in os.listdir(path) if tag in f]
    mjd = []
    for file in files:
        new_path = os.path.join(path, file)
        mjd.append(fits.getheader(new_path)['MJD'])
    
    return [x for _, x in sorted(zip(mjd, files))]

def _find_img_closest_value(path:str, files:list, keyword:str, value):
    mjds = []
    for file in files:
        file = os.path.join(path, file)
        mjd = fits.getheader(file)[keyword]
        mjds.append(mjd)

    index = np.argmin(np.abs(np.array(mjds)-value))
    return files[index]

def get_coords_in_series(path:str, dates:list, mjds:list):
    xcoord, ycoord, header_dates = np.array([]), np.array([]), []
    files = os.listdir(path)
    for file in files:
        hdr = fits.getheader(os.path.join(path, file))
        header_dates.append(hdr['date'])

    for date , mjd in zip(dates, mjds):
        indexes = np.where(np.array(header_dates) == date)[0]
        new_files = [files[idx] for idx in indexes]
        file = _find_img_closest_value(path, new_files, 'MJD', mjd)
        x, y, _ = get_obj_coords(path, files[indexes[-1]])
        xcoord = np.append(xcoord, x)
        ycoord = np.append(ycoord, y)
    return xcoord, ycoord