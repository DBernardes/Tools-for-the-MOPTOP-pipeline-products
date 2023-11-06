# # -*- coding: utf-8 -*-

import pytest
import os
import astropy.io.fits as fits
import pandas as pd
import numpy as np
from photometry import Photometry, Object
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from scipy.interpolate import UnivariateSpline
from copy import copy


@pytest.fixture()
def obj():
    return Object(
        "target",
        xcoord=100,
        ycoord=200,
        mjd=300,
        ra="00:00:00",
        dec="00:00:00",
        sky_photons=15,
        star_photons=1000,
        star_err=10,
        star_mag=14,
        star_mag_err=0.1,
    )


def test_object_initialization(obj):
    assert obj.name == "target"
    assert obj.xcoord == 100
    assert obj.ycoord == 200
    assert obj.mjd == 300
    assert obj.ra == "00:00:00"
    assert obj.dec == "00:00:00"
    assert obj.sky_photons == 15
    assert obj.star_photons == 1000
    assert obj.star_err == 10
    assert obj.star_mag == 14
    assert obj.star_mag_err == 0.1


# ----------------------------------------------------------------

base_path = os.path.join("tests", "images")
file_path = os.path.join(base_path, "3_e_20230818_5_1_1_1.fits")
obj_coordinates = os.path.join(base_path, "..", "setup", "objects coordinates.csv")
max_size = 10
bkg_sigma = 4
star_name = "comparison"
max_radius = max_size // 2
image, header = fits.getdata(file_path, header=True)
image *= header["GAIN"]
image_shape = image.shape
ra, dec = "19:03:40.0436", "+40:50:28.191"
xcoord, ycoord = 736, 624
obj_list = [Object(star_name, xcoord, ycoord, header["mjd"])]


coordinates = {"name": [star_name], "ra_cam3": [ra], "dec_cam3": [dec]}
coordinates = pd.DataFrame.from_dict(coordinates)


@pytest.fixture()
def phot():
    return Photometry(file_path, coordinates, max_size)


def test_convert_world_to_pixel(phot):
    tup = phot._convert_coords_to_pixel(ra, dec)

    assert (xcoord, ycoord) == tup


def test_initialization(phot):
    assert phot.fits_file == file_path
    assert phot.max_radius == max_radius
    assert phot.image.all() == image.all()
    assert phot.image_shape == image_shape
    assert phot.header == header
    assert phot.obj_list == obj_list


def test_convert_pixels_to_world(phot):
    tup = phot._convert_pixels_to_coords(xcoord, ycoord)

    wcs = WCS(header)
    ra, dec = wcs.pixel_to_world(xcoord, ycoord).to_string("hmsdms").split(" ")
    assert tup == (ra, dec)


size = max_radius
new_image = copy(image[ycoord - size : ycoord + size, xcoord - size : xcoord + size])


def test_calculate_bkg_level(phot):
    background_level = phot._calc_background_level(new_image, bkg_sigma)

    median = np.median(new_image)
    std = np.median(np.abs(new_image - median))
    max_lim, min_lim = median + bkg_sigma * std, median - bkg_sigma * std
    indexes = np.where((min_lim < new_image) & (new_image < max_lim))
    background_pixels = new_image[indexes]
    sky = np.median(background_pixels)

    assert background_level == sky + 4 * np.std(background_pixels)


def test_create_obj_mask(phot):
    phot.bkg_sigma = 4
    new_working_mask = phot._create_objects_mask(image)

    background_level = phot._calc_background_level(image, bkg_sigma)
    y_coords, x_coords = np.where(image > background_level)
    working_mask = np.zeros(image.shape, bool)
    working_mask[(y_coords, x_coords)] = 1

    for x, y in zip(x_coords, y_coords):
        sub_image = working_mask[y - 1 : y + 2, x - 1 : x + 2]
        if np.sum(sub_image) <= 1:
            working_mask[y, x] = 0

    assert working_mask.all() == new_working_mask.all()


def test_find_closest_bject(phot):
    phot.bkg_sigma = 4
    tup = phot._find_closest_object(new_image)

    wk_mask = phot._create_objects_mask(new_image)
    center_coord = wk_mask.shape[0] // 2
    y_coords, x_coords = np.where(wk_mask == 1)
    x_coords -= center_coord
    y_coords -= center_coord
    distances = np.sqrt(x_coords**2 + y_coords**2)
    idx_min = np.argmin(distances)
    closest_x = x_coords[idx_min] + center_coord
    closest_y = y_coords[idx_min] + center_coord

    assert tup == (closest_x, closest_y)


def test_find_coords_max_pixel(phot):
    phot.bkg_sigma = 4
    tup = phot._find_coords_max_pixel(new_image)

    working_mask = phot._create_objects_mask(new_image)
    max_value = np.max(new_image)
    if np.sum(working_mask) > 0:
        max_value = np.max(new_image[working_mask])
    new_y, new_x = np.where(new_image == max_value)

    assert tup == (new_x[0], new_y[0])


def test_reset_obj_coords(phot):
    phot.reset_object_coords()
    obj = phot.obj_list[0]
    assert obj.xcoord == 744
    assert obj.ycoord == 625


def test_calc_estimate_sky_photons(phot):
    new_sky = phot._calc_estimate_sky_photons(new_image)

    median = np.median(new_image)
    std = np.median(np.abs(new_image - median))

    max_lim, min_lim = median + bkg_sigma * std, median - bkg_sigma * std
    indexes = np.where((min_lim < new_image) & (new_image < max_lim))
    background_pixels = new_image[indexes]
    sky = np.median(background_pixels)

    assert sky == new_sky


def test_calculate_star_radius(phot):
    phot.calculate_star_radius()

    _object = obj_list[0]
    x, y = _object.xcoord, _object.ycoord
    r = max_radius
    img_data = copy(image[y - r : y + r, x - r : x + r])
    sky_photons = phot._calc_estimate_sky_photons(img_data)
    img_data -= sky_photons

    light_profile = np.take(img_data, r - 1, axis=0)
    half_max = np.max(light_profile) / 2
    n = len(light_profile)
    x = np.linspace(0, n - 1, n)
    spline = UnivariateSpline(x, light_profile - half_max, s=None)

    roots = spline.roots()
    idx_max_val = np.argmax(spline(x))
    tmp = np.abs(roots - idx_max_val)
    idx_r_1 = np.argmin(tmp)
    tmp[idx_r_1] = 1e10
    idx_r_2 = np.argmin(tmp)

    fwhm = np.abs(roots[idx_r_2] - roots[idx_r_1])
    assert phot.star_radius == 3 * fwhm


def test_create_sky_mask(phot):
    phot.reset_object_coords()
    phot.calculate_star_radius()
    obj = phot.obj_list[0]
    star_radius = phot.star_radius
    mask = phot._create_sky_mask(obj.xcoord, obj.ycoord, star_radius)

    working_mask = np.ones(image_shape, bool)
    ym, xm = np.indices(image_shape, dtype="float32")
    r = np.sqrt((xm - xcoord) ** 2 + (ym - ycoord) ** 2)
    new_mask = (r > 2 * star_radius) * (r < 3 * star_radius) * working_mask

    assert mask.all() == new_mask.all()


def test_calculate_sky_photons(phot):
    phot.reset_object_coords()
    obj = phot.obj_list[0]
    phot.calculate_star_radius()
    star_radius = phot.star_radius
    phot.calc_sky_photons()

    mask = phot._create_sky_mask(obj.xcoord, obj.ycoord, star_radius)
    sky = image[np.where(mask)]
    sky_photons = np.median(sky)

    assert sky_photons == phot.obj_list[0].sky_photons


def test_creat_star_mask(phot):
    phot.reset_object_coords()
    phot.calculate_star_radius()
    star_radius = phot.star_radius
    obj = phot.obj_list[0]
    mask = phot._create_star_mask(obj.xcoord, obj.ycoord, star_radius)

    working_mask = np.ones(image_shape, bool)
    ym, xm = np.indices(image_shape, dtype="float32")
    r = np.sqrt((xm - xcoord) ** 2 + (ym - ycoord) ** 2)
    new_mask = (r < star_radius) * working_mask

    assert mask.all() == new_mask.all()


def test_calc_star_photons(phot):
    phot.reset_object_coords()
    phot.calculate_star_radius()
    phot.calc_sky_photons()
    phot.calc_star_photons()
    star_photons = phot.obj_list[0].star_photons

    obj = phot.obj_list[0]
    star_radius = phot.star_radius
    sky_photons = obj.sky_photons

    mask = phot._create_star_mask(obj.xcoord, obj.ycoord, star_radius)
    star = image[np.where(mask)]
    new_star_photons = np.sum(star - sky_photons)

    star_err = np.sqrt(star_photons + sky_photons * star.shape[0])

    assert new_star_photons == obj.star_photons
    assert star_err == obj.star_err
