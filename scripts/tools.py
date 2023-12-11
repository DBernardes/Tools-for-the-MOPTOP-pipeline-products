#!/usr/bin/env python

"""tools.py: this file has a set of functions developed for the analysis of the products obtained running the MOPTOP pipeline in a data set."""

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import pandas as pd
import numpy as np
import astropy.io.fits as fits
import shutil
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from numpy import ndarray

import matplotlib.pyplot as plt

# (pol-%, phase-deg)
low_polarized_stars = {
    "HD14069": {"B": (0.111, 93.62), "V": (0.022, 156.57)},
    "BD+32 3739": {"B": (0.039, 79.38), "V": (0.025, 35.79)},
}
# https://www.not.iac.es/instruments/turpol/std/zpstd.html
high_polarized_stars = {
    "BD+64 106": 0,
    "HD251204": 0,
    "VICyg12": 0,
    "HD155197": 0,
    "HILT960": 0,
    "GRB 230818A": 0,
}
# https://www.not.iac.es/instruments/turpol/std/hpstd.html
LIMIT_MJD = 59774

parameters_MAS_estimator = {
    "p_min": (1, 0.72, 0.6, -0.83, 4.41),
    "p_max": (1, 0.97, 2.01),
}


def manipulate_csv_file(
    path: str, dest_path: str = "", filename: str = "manipulated_data"
):
    """Manipulate the raw_data.csv file created by the MOPTOP pipeline

    Parameters
    ----------
    path : str
        Path of the folder where the raw_data.csv file is
    dest_path : str
        Path of the folder where the raw_data.csv file will be saved
    filename : str
        Name of the new .csv file
    """
    df = pd.read_csv(path)
    df = df[df["q_avg"].notna()]
    for column_name in df.columns.values:
        if column_name not in [
            "q_avg",
            "q_err",
            "u_avg",
            "u_err",
            "date",
            "mjd",
            "wave",
            "s1_src_cam1",
            "s1_src_cam1_err",
            "s1_src_cam2",
            "s1_src_cam2_err",
            "time",
            "x",
            "y",
        ]:
            df = df.drop(column_name, axis=1)
    df = df.sort_values(by=["wave", "mjd"])
    if dest_path == "":
        base_path = path.split("/")[:-1]
        dest_path = os.path.join(*base_path)
    file_path = os.path.join(dest_path, filename + ".csv")
    pd.DataFrame.to_csv(df, file_path, index=False)


def calculate_polarization_1(path: str) -> dict:
    """Calculate the polarization for a set of q and u values

    Parameters
    ----------
    path : str
        path of the manipulated_data.csv file

    Returns
    -------
    dict
        polarization values
    """
    qu_dict = sort_qu_per_filter(path)
    pol_dict = {"B": [], "V": [], "R": [], "I": [], "L": []}
    for key in qu_dict.keys():
        mjd = np.asarray(qu_dict[key]["mjd"])
        q = np.asarray(qu_dict[key]["q"])
        u = np.asarray(qu_dict[key]["u"])
        q_err = np.asarray(qu_dict[key]["q_err"])
        u_err = np.asarray(qu_dict[key]["u_err"])
        pol = np.sqrt(q**2 + u**2) * 100
        err = np.sqrt((q / pol) ** 2 * q_err**2 + (u / pol) ** 2 * u_err**2) * 100
        pol_dict[key].append(mjd)
        pol_dict[key].append(pol)
        pol_dict[key].append(err)

    return pol_dict


def calculate_polarization_and_phase(
    q: ndarray, q_err: ndarray, u: ndarray, u_err: ndarray
) -> tuple[ndarray, ndarray]:
    """Calculates the polarization and the phase for a set of q and u values

    Parameters
    ----------
    q : ndarray
        q Stokes values
    q_err : ndarray
        error in q Stokes values
    u : ndarray
        u Stokes values
    u_err : ndarray
        u error in u Stokes values

    Returns
    -------
    tuple[ndarray, ndarray]
        polarization and phase values, with their respective errors
    """
    pol = np.sqrt(q**2 + u**2) * 100
    pol_err = np.sqrt((q / pol) ** 2 * q_err**2 + (u / pol) ** 2 * u_err**2) * 100

    # * Propagating the errors
    x = u / q
    std_x = np.abs(x) * np.sqrt((q_err / q) ** 2 + (u_err / u) ** 2)

    phase = np.rad2deg(0.5 * np.arctan(x))
    phase_err = np.rad2deg(0.5 * std_x / (1 + x**2))

    return pol, pol_err, phase, phase_err


def sort_qu_per_filter(path: str) -> dict:
    """Sort the q and u values for the filter BVRIL

    Parameters
    ----------
    path : str
        path for the manipulated_data.csv file

    Returns
    -------
    dict
        q and u values sorted for the filters BVRIL
    """
    df = pd.read_csv(path)
    qu_dict = {}
    for filter in ["B", "V", "R", "I", "L"]:
        rows = df.loc[df["wave"] == f"MOP-{filter}"]
        qu_dict[filter] = {
            "mjd": rows["mjd"],
            "q": rows["q_avg"],
            "u": rows["u_avg"],
            "q_err": rows["q_err"],
            "u_err": rows["u_err"],
        }

    return qu_dict


def track_obj_over_images(
    path: str, tag: str = ".fits"
) -> tuple[ndarray, ndarray, ndarray]:
    """Track an object over a series of images

    Parameters
    ----------
    path : str
        path of the folder with the images
    tag : str, optional
        a tag to be used for the image names when listing the folder content, by default ".fits"

    Returns
    -------
    tuple[ndarray, ndarray, ndarray]
        Arrays of the x and y coordinates of the object, as well as the MJD values obtained for the set of images
    """
    files = sort_files(path, tag)
    xcoord, ycoord, mjds = np.array([]), np.array([]), np.array([])
    for file in files:
        file = os.path.join(path, file)
        x, y, mjd = get_obj_coords(file)
        xcoord = np.append(xcoord, x)
        ycoord = np.append(ycoord, y)
        mjds = np.append(mjds, mjd)
    return xcoord, ycoord, mjds


def get_obj_coords(
    file: str, ra: str = None, dec: str = None, unit=(u.hourangle, u.deg)
) -> tuple[int, int, float]:
    """Get the object coordinates in image.

    Parameters
    ----------
    file : str
        file name
    ra : str, optional
        object right ascension, by default None.
    dec : str, optional
        object declination, by default None.

    Returns
    -------
    tuple[int, int, float]
        X and Y coordinates of the object, together with the MJD of the image.
    """
    hdr = fits.getheader(file)
    wcs = WCS(hdr)
    coords_str = f'{hdr["CAT-RA"]} {hdr["CAT-DEC"]}'
    if None not in [ra, dec]:
        coords_str = f"{ra} {dec}"
    coord = SkyCoord(coords_str, frame="fk5", unit=unit)
    x, y = wcs.world_to_pixel(coord)
    return int(x) + 1, int(y) + 1, hdr["MJD"]


def select_images_keyword_interval(
    path: str,
    dest_path: str,
    keyword: str,
    min: float = -np.infty,
    max: float = np.infty,
    tag: str = ".fits",
) -> None:
    """Select images based on an interval of allowed values for a keyword

    Parameters
    ----------
    path : str
        path of the folder with the set of images
    dest_path : str
        destination path
    keyword : str
        name of the header keyword to be used
    min : float, optional
        minimum value of the keyword, by default -np.infty
    max : float, optional
        maximum value of the keyword, by default np.infty
    tag : str, optional
        a tag to be used for the image name when listing the folder content, by default ".fits"
    """
    files = sort_files(path, tag)
    os.makedirs(dest_path, exist_ok=True)
    for file in files:
        src_file = os.path.join(path, file)
        hdr = fits.getheader(src_file)
        if min < hdr[keyword] < max:
            dest_file = os.path.join(dest_path, file)
            shutil.copyfile(src_file, dest_file)


def select_images_keyword_value(
    path: str,
    dest_path: str,
    keyword: str,
    value: float | int | str,
    tag: str = ".fits",
    sort_images=False,
) -> None:
    """Select images based on the value of a keyword

    Parameters
    ----------
    path : str
        folder with the set of images
    dest_path : str
        destination path
    keyword : str
        keyword name
    value : float | int | str
        keyword value
    tag : str, optional
        a tag to be used for the image names when listing the folder content, by default ".fits"
    sort_images : bool, optional
        if the function should or not to sort the images based on their MHD.
        This procedure could take a little more time to be done. By default False
    """
    if sort_images:
        files = sort_files(path, tag)
    else:
        files = [f for f in os.listdir(path) if tag in f]
    os.makedirs(dest_path, exist_ok=True)
    for file in files:
        print(file)
        src_file = os.path.join(path, file)
        hdr = fits.getheader(src_file)
        if hdr[keyword] == value:
            dest_file = os.path.join(dest_path, file)
            shutil.copyfile(src_file, dest_file)


def delete_file_keyword_value(
    path: str, keyword: str, value: float | str | int, tag: str = ".fits"
) -> None:
    """Delete images based on a keyword

    Parameters
    ----------
    path : str
        path of the folder with the set of images
    keyword : str
        keyword name
    value : float | str | int
        keyword value
    tag : str, optional
        a tag to be used for the image names when listing the folder content, by default ".fits"
    """
    files = sort_files(path, tag)
    for file in files:
        src_file = os.path.join(path, file)
        hdr = fits.getheader(src_file)
        if hdr[keyword] == value:
            os.remove(src_file)
    return


def calculate_mean_images(path: str, tag: str = ".fits") -> ndarray:
    """Calculate the mean value of a set of images

    Parameters
    ----------
    path : str
        path of the folder with the set of images
    tag : str, optional
        a tag to be used for the image names when listing the folder content, by default ".fits"

    Returns
    -------
    ndarray
        mean value calculate for the set of images
    """
    files = sort_files(path, tag)
    mean = []
    for file in files:
        image = fits.getdata(os.path.join(path, file))
        tmp = np.mean(image)
        mean.append(tmp)
    return np.asarray(mean)


def calculate_maximum_images(path, tag: str = ".fits") -> ndarray:
    """Calculate the maximum value of a set of images

    Parameters
    ----------
    path : str
        path of the folder with the set of images
    tag : str, optional
        a tag to be used for the image names when listing the folder content, by default ".fits"

    Returns
    -------
    ndarray
        maximum value calculate for the set of images
    """
    files = sort_files(path, tag)
    _max = []
    for file in files:
        image = fits.getdata(os.path.join(path, file))
        tmp = np.max(image)
        _max.append(tmp)
    return np.asarray(_max)


def sort_files(path: str, tag: str = ".fits"):
    files = [f for f in os.listdir(path) if tag in f]
    mjd = []
    for file in files:
        new_path = os.path.join(path, file)
        mjd.append(fits.getheader(new_path)["MJD"])

    return [x for _, x in sorted(zip(mjd, files))]


def _find_img_closest_value(path: str, files: list, keyword: str, value):
    mjds = []
    for file in files:
        file = os.path.join(path, file)
        mjd = fits.getheader(file)[keyword]
        mjds.append(mjd)

    index = np.argmin(np.abs(np.array(mjds) - value))
    return files[index]


def get_coords_in_series(
    path: str, dates: list, mjds: list, tag: str = ".fits"
) -> tuple[float, float]:
    """Get the x and y coordinates of the object in a series of images found in path

    Parameters
    ----------
    path : str
        path of the folder with the set of images
    dates : list
        dates of those images in which the coordinates should be calculated
    mjds : list
        mjd of those images in which the coordinates should be calculated
    tag : str, optional
        a tag for the files in the folder that should be used, by default ".fits"

    Returns
    -------
    tuple[float, float]
        _description_
    """
    xcoord, ycoord, header_dates = np.array([]), np.array([]), []
    files = [f for f in os.listdir(path) if tag in f]
    for file in files:
        hdr = fits.getheader(os.path.join(path, file))
        header_dates.append(hdr["date"])

    for date, mjd in zip(dates, mjds):
        indexes = np.where(np.array(header_dates) == date)[0]
        new_files = [files[idx] for idx in indexes]
        file = _find_img_closest_value(path, new_files, "MJD", mjd)
        x, y, _ = get_obj_coords(path, file)
        xcoord = np.append(xcoord, x)
        ycoord = np.append(ycoord, y)
    return xcoord, ycoord


def sigma_clipping(x: list, y: list, x_err=[], y_err=[], sigma=5, iter=1) -> tuple:
    """Perform a sigma clipping for the x and y parameters

    Parameters
    ----------
    x : list
        any parameter to be clipped
    y : list
        any parameter to be clipped
    x_err : list, optional
        standard deviation of the x paramrter, by default []
    y_err : list, optional
        standard deviation of the y parameter, by default []
    sigma : int, optional
        number of standard deviations to be used for clipping, by default 5
    iter : int, optional
        number of iterations to be used in clipping, by default 1

    Returns
    -------
    tuple
        clipped x and y parameters, with their corresponding standard deviation
    """
    (
        x,
        y,
    ) = np.asarray(
        x
    ), np.asarray(y)
    if len(x_err) == 0:
        x_err = np.zeros(len(x))
    else:
        x_err = np.asarray(x_err)
    if len(y_err) == 0:
        y_err = np.zeros(len(y))
    else:
        y_err = np.asarray(y_err)
    for _ in range(iter):
        medianx, mediany = np.median(x), np.median(y)
        stdx = np.median(np.abs(x - medianx))
        stdy = np.median(np.abs(y - mediany))
        indexes = np.where((medianx - sigma * stdx < x) & (x < medianx + sigma * stdx))
        x, y = x[indexes], y[indexes]
        x_err, y_err = x_err[indexes], y_err[indexes]
        indexes = np.where((mediany - sigma * stdy < y) & (y < mediany + sigma * stdy))
        x, y = x[indexes], y[indexes]
        x_err, y_err = x_err[indexes], y_err[indexes]

    return x, y, x_err, y_err


def read_calculated_qu_values() -> dict:
    """Read the calculated q and u values found in the 'mean_qu_values' csv file.

    Returns
    -------
    dict
        calculated q and u values for the filter UBVRI of MOPTOP.
    """
    caculated_qu = {
        "B": (),
        "V": (),
        "R": (),
        "I": (),
        "L": (),
    }
    csv_file_name = os.path.join("csv", "mean_qu_values.csv")
    df = pd.read_csv(csv_file_name)
    for star in df["star"]:
        if star not in ["GD319", "HD14069", "BD+32 3739"]:
            df.drop(df[df["star"] == star].index, inplace=True)
    for filter in caculated_qu.keys():
        rows = df.loc[df["filter"] == filter]
        q, u = rows["q"], rows["u"]
        caculated_qu[filter] = (np.mean(q), np.std(q), np.mean(u), np.std(u))
    return caculated_qu


def get_instrumental_polarization(
    x: ndarray, y: ndarray, _filter: str, parameter: str
) -> ndarray:
    """Get the instrumental polarization for a given filter and Stokes parameter

    Parameters
    ----------
    x : ndarray
        X coordinates of the object
    y : ndarray
        Y coordinates of the object
    _filter : str
        Johnson filter used in observation
    parameter : str
        q or u Stokes parameter

    Returns
    -------
    ndarray
        q or u instrumental polarization values
    """
    csv_path = os.path.join(
        "scripts",
        "csv",
        "plane_coefficients.csv",
    )
    df = pd.read_csv(csv_path)
    row = df.loc[df["filter"] == _filter]

    a, b, c = [row[f"{parameter}{letter}"].values[0] for letter in ["a", "b", "c"]]
    Z = a * x + b * y + c
    return Z


def novel_pol_error(
    q: ndarray, q_err: ndarray, u: ndarray, u_err: ndarray
) -> tuple[ndarray]:
    """Calculate the novel polarization

    Parameters
    ----------
    q : ndarray
        q Stokes parameter
    q_err : ndarray
        q error
    u : ndarray
        u Stokes parameter
    u_err : ndarray
        u error

    Returns
    -------
    tuple[ndarray]
        The values of the naive, modified, minimum, and maximum polarizations, and standard deviations.
    """
    pol = np.sqrt(q**2 + u**2) * 100  # eq2
    pol_err = np.sqrt((q_err**2 + u_err**2) / 2) * 100  # eq29

    b_sq = (q**2 * u_err**2 + u**2 * q_err**2) / (q**2 + u**2)  # eq30
    p_mas_err = (
        np.sqrt((u**2 * u_err**2 + q**2 * q_err**2) / (q**2 + u**2)) * 100
    )  # eq31
    p_mas = pol - (b_sq / (2.0 * pol)) * (1.0 - np.exp(-(pol**2) / b_sq))  # eq37

    # eq26
    p_alpha, beta, gamma = parameters_MAS_estimator["p_max"]
    p_max = p_mas + pol_err * p_alpha * (1.0 - beta * np.exp(-gamma * p_mas / pol_err))
    p_alpha, beta, gamma, omega, phi = parameters_MAS_estimator["p_min"]
    p_min = p_mas - pol_err * p_alpha * (
        1.0
        + beta
        * np.exp(-gamma * p_mas / pol_err)
        * np.sin(omega * p_mas / pol_err + phi)
    )

    return (
        pol,
        pol_err,
        p_mas,
        p_mas_err,
        p_min,
        p_max,
    )


def calculate_qu_low_standards(star: str, filter: str):
    if filter in ["R", "I", "L"]:
        filter = "V"
    pol, phase = low_polarized_stars[star][filter]
    pol, phase = pol / 100, np.deg2rad(phase)
    q = pol / np.sqrt(1 + np.tan(2 * phase) ** 2)
    u = np.sqrt(pol**2 - q**2)
    return (q, u)
