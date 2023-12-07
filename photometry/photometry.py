#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"

import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import UnivariateSpline
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from dataclasses import dataclass
from pandas import DataFrame
import matplotlib.pyplot as plt
from math import sqrt, log10, log2
from sbpy.calib import vega_fluxd
from scipy.constants import h, c
from copy import copy


class Photometry:
    """Photometry class"""

    READ_NOISE = 0.8  # e-
    DARK_CURRENT = 0.3  # e-/pix/s

    def __init__(self, fits_file: str, objects: DataFrame, max_size: int = 15) -> None:
        """Initialize the class

        Parameters
        ----------
        fits_file : str
            FITS file name.
        objects: DataFrame
            A pandas Dataframe of tuples with the name, right ascension, and declination of the objects.
        max_size: int
            maximum size, in pixels, of the box in which the object can be found. Default to 30.
        """
        self.fits_file = fits_file
        self.max_radius = max_size // 2
        self.image, self.header = fits.getdata(fits_file, header=True)
        self.image *= self.header["GAIN"]
        self.image_shape = self.image.shape
        self.obj_list = []
        for _object in objects.itertuples(name=None, index=False):
            name, ra, dec = _object
            xcoord, ycoord = self._convert_coords_to_pixel(ra, dec)
            self.obj_list.append(Object(name, xcoord, ycoord, self.header["mjd"]))
        return

    def _convert_coords_to_pixel(
        self,
        ra: str,
        dec: str,
    ) -> tuple[int, int]:
        """Get the object coordinates in image.

        Parameters
        ----------
        ra : str, optional
            object right ascension, by default None.
        dec : str, optional
            object declination, by default None.

        Returns
        -------
        tuple[int, int, float]
            X and Y coordinates of the object, together with the MJD of the image.
        """

        wcs = WCS(self.header)
        coords_str = f"{ra} {dec}"
        coord = SkyCoord(coords_str, frame="fk5", unit=(u.hourangle, u.deg))
        x, y = wcs.world_to_pixel(coord)
        if x < 0 or y < 0:
            raise ValueError(f"Object is out of field: ({x},{y}).")
        xcoord = int(x) + 1
        ycoord = int(y) + 1
        return xcoord, ycoord

    def _convert_pixels_to_coords(self, x, y) -> tuple[str, str]:
        wcs = WCS(self.header)

        ra, dec = (
            wcs.pixel_to_world(
                x,
                y,
            )
            .to_string("hmsdms")
            .split(" ")
        )
        return ra, dec

    @staticmethod
    def _calc_background_level(image: np.ndarray, bkg_sigma: float) -> float:
        median = np.median(image)
        std = np.median(np.abs(image - median))

        max_lim, min_lim = median + bkg_sigma * std, median - bkg_sigma * std
        indexes = np.where((min_lim < image) & (image < max_lim))
        background_pixels = image[indexes]
        sky = np.median(background_pixels)
        background_level = sky + 4 * np.std(background_pixels)

        return background_level

    def _create_objects_mask(self, image: np.ndarray) -> np.ndarray:
        background_level = self._calc_background_level(image, self.bkg_sigma)

        working_mask = np.zeros(image.shape, bool)
        y_coords, x_coords = np.where(image > background_level)
        working_mask[(y_coords, x_coords)] = 1

        for x, y in zip(x_coords, y_coords):
            sub_image = working_mask[y - 1 : y + 2, x - 1 : x + 2]
            if np.sum(sub_image) <= 1:
                working_mask[y, x] = 0

        return working_mask

    def _find_closest_object(self, image):
        working_mask = self._create_objects_mask(image)
        if np.sum(working_mask) > 0:
            center_coord = working_mask.shape[0] // 2
            y_coords, x_coords = np.where(working_mask == 1)
            x_coords -= center_coord
            y_coords -= center_coord
            distances = np.sqrt(x_coords**2 + y_coords**2)
            idx_min = np.argmin(distances)
            closest_x = x_coords[idx_min] + center_coord
            closest_y = y_coords[idx_min] + center_coord

            return closest_x, closest_y
        else:
            return None

    def _find_coords_max_pixel(self, image: np.ndarray):
        working_mask = self._create_objects_mask(image)

        max_value = np.max(image)
        if np.sum(working_mask) > 0:
            max_value = np.max(image[working_mask])
        new_y, new_x = np.where(image == max_value)

        return new_x[0], new_y[0]

    def reset_object_coords(self, bkg_sigma: int = 4):
        """Recalculate the object coordinates."""
        size = self.max_radius
        self.bkg_sigma = bkg_sigma
        for idx, _object in enumerate(self.obj_list):
            x, y = _object.xcoord, _object.ycoord
            image = self.image[y - size : y + size, x - size : x + size]

            new_coords = self._find_closest_object(image)
            if new_coords != None:
                closest_x, closest_y = new_coords
                x += closest_x - size
                y += closest_y - size
                image = self.image[
                    y - size : y + size,
                    x - size : x + size,
                ]

            new_x, new_y = self._find_coords_max_pixel(image)
            new_x += x - size
            new_y += y - size
            ra, dec = self._convert_pixels_to_coords(new_x, new_y)
            _object.ra, _object.dec = ra, dec
            _object.xcoord, _object.ycoord = new_x + 1, new_y + 1

            self.obj_list[idx] = _object
        return

    def _calc_estimate_sky_photons(self, image, bkg_sigma: float = 4):
        # This is a first estimate for the number of photons of the sky

        median = np.median(image)
        std = np.median(np.abs(image - median))

        max_lim, min_lim = median + bkg_sigma * std, median - bkg_sigma * std
        indexes = np.where((min_lim < image) & (image < max_lim))
        background_pixels = image[indexes]
        sky = np.median(background_pixels)

        return sky

    def calculate_star_radius(self, coeff_radius_fwhm: float = 3) -> float:
        """Calculate FWHM of the object.

        Parameters
        ----------
        coeff_radius_fwhm: float
            The multiplicative coefficient between the star radius and the FWHM.

        Returns
        -------
        star_radius: float
            The calculated FWHM for the object, times the coefficient_radius_fwhm parameter.
        """

        _object = [obj for obj in self.obj_list if "candidate1" in obj.name][0]
        x, y = _object.xcoord, _object.ycoord
        r = self.max_radius
        img_data = copy(self.image[y - r : y + r, x - r : x + r])
        # plt.imshow(img_data, origin="lower")
        # plt.show()
        sky_photons = self._calc_estimate_sky_photons(img_data)
        img_data -= sky_photons

        light_profile = np.take(img_data, r - 1, axis=0)
        # plt.plot(light_profile)
        # plt.show()
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

        # plt.plot(x, spline(x), "b-")
        # plt.plot(roots[idx_r_2], spline(roots[idx_r_2]), "bo")
        # plt.plot(roots[idx_r_1], spline(roots[idx_r_1]), "bo")
        # print(roots)
        # plt.show()

        fwhm = np.abs(roots[idx_r_2] - roots[idx_r_1])
        self.star_radius = coeff_radius_fwhm * fwhm

        return self.star_radius

    def _create_sky_mask(self, xcoord, ycoord, star_radius):
        working_mask = np.ones(self.image_shape, bool)
        ym, xm = np.indices(self.image_shape, dtype="float32")
        r = np.sqrt((xm - xcoord) ** 2 + (ym - ycoord) ** 2)
        mask = (r > 2 * star_radius) * (r < 3 * star_radius) * working_mask
        return mask

    def calc_sky_photons(self):
        """Calculate the number of photons of the sky"""
        for idx, _object in enumerate(self.obj_list):
            mask = self._create_sky_mask(
                _object.xcoord, _object.ycoord, self.star_radius
            )
            sky = self.image[np.where(mask)]
            self.obj_list[idx].sky_photons = np.median(sky)

    def _create_star_mask(self, xcoord, ycoord, star_radius):
        working_mask = np.ones(self.image_shape, bool)
        ym, xm = np.indices(self.image_shape, dtype="float32")
        r = np.sqrt((xm - xcoord) ** 2 + (ym - ycoord) ** 2)
        mask = (r < star_radius) * working_mask
        return mask

    def calc_star_photons(self):
        """Calculate the number of photons of the object."""
        t_exp = self.header["EXPTIME"]
        for idx, _object in enumerate(self.obj_list):
            mask = self._create_star_mask(
                _object.xcoord, _object.ycoord, self.star_radius
            )
            star = self.image[np.where(mask)]
            star_photons = np.sum(star - _object.sky_photons)
            if star_photons < 0:
                star_photons = 0
            self.obj_list[idx].star_photons = star_photons
            self.obj_list[idx].star_err = np.sqrt(
                star_photons
                + star.shape[0]
                * (
                    _object.sky_photons
                    + t_exp * self.DARK_CURRENT
                    + self.READ_NOISE**2
                )
            )
        return

    def calc_magnitude(self):
        """Calculate the magnitude for each object in the list"""
        for obj in self.obj_list:
            _filter, exptime = self.header["FILTER1"][-1], self.header["EXPTIME"]
            obj.calc_magnitude(_filter, exptime)


@dataclass
class Object:
    """Class to keep the photometry information related to an astronomic object."""

    name: str
    xcoord: str
    ycoord: str
    mjd: float
    ra: str = 0
    dec: str = 0
    sky_photons: float = 0
    star_photons: float = 0
    star_err: float = 0
    star_mag: float = 0
    star_mag_err: float = 0

    _VEGA = vega_fluxd.get()
    _LT_EFFECTIVE_AREA = 2.982  # m2

    def calc_magnitude(self, _filter: str, exp_time: float):
        """Calculate the object magnitude

        Parameters
        ----------
        _filter : str
            UBVRI filter used in observation.
        exp_time : float
            Exposure time
        """

        keyword_str = "Johnson"
        if _filter not in ["U", "B", "V"]:
            keyword_str = "Cousins"

        eff_lambda = self._VEGA[f"{keyword_str} {_filter}(lambda eff)"].value * 1e-6
        photon_energy = h * c / eff_lambda
        vega_photons = (
            self._VEGA[f"{keyword_str} {_filter}"].value
            * 1e7
            * self._LT_EFFECTIVE_AREA
            * eff_lambda
            * exp_time
            / photon_energy
        )

        self.star_mag_err = 2.5 * self.star_err / (self.star_photons * log2(10))
        self.star_mag = -2.5 * log10(self.star_photons / vega_photons)
