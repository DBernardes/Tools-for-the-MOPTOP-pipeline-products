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


class Photometry:
    """Photometry class"""

    def __init__(self, file: str, objects: DataFrame, max_size: int = 15) -> None:
        """Initialize the class

        Parameters
        ----------
        file : str
            FITS file name.
        objects: DataFrame
            A pandas Dataframe of tuples with the name, right ascension, and declination of the objects.
        max_size: int
            maximum size, in pixels, of the box in which the object can be found. Default to 30.
        """
        self.file = file
        self.max_radius = max_size // 2
        self.image, self.header = fits.getdata(file, header=True)
        self.image *= self.header["GAIN"]  # ? is this right
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
    def _calc_background_level(image: np.ndarray, nsigma: int = 4) -> float:
        median = np.median(image)
        std = np.median(np.abs(image - median))

        indexes = np.where(
            (median - nsigma * std < image) & (image < median + nsigma * std)
        )
        background_pixels = image[indexes]
        sky = np.median(background_pixels)
        sky_err = np.std(background_pixels)

        return sky, sky_err

    @staticmethod
    def _create_objects_mask(image, background_level: float) -> np.ndarray:
        working_mask = np.zeros(image.shape, bool)
        y_coords, x_coords = np.where(image > background_level)
        working_mask[(y_coords, x_coords)] = 1

        for x, y in zip(x_coords, y_coords):
            sub_image = working_mask[y - 1 : y + 2, x - 1 : x + 2]
            if np.sum(sub_image) <= 1:
                working_mask[y, x] = 0

        # plt.imshow(working_mask, origin="lower")
        # plt.show()

        return working_mask

    @staticmethod
    def _find_closest_object(working_mask: np.ndarray):
        center_coord = working_mask.shape[0] // 2
        y_coords, x_coords = np.where(working_mask == 1)
        x_coords -= center_coord
        y_coords -= center_coord
        distances = np.sqrt(x_coords**2 + y_coords**2)
        idx_min = np.argmin(distances)
        closest_x = x_coords[idx_min] + center_coord
        closest_y = y_coords[idx_min] + center_coord

        return closest_x, closest_y

    def _find_coords_max_pixel(self, closest_x, closest_y):
        size = 10

        image = self.image[
            closest_y - size : closest_y + size, closest_x - size : closest_x + size
        ]

        sky, sky_err = self._calc_background_level(image)
        background_level = sky + 4 * sky_err
        working_mask = self._create_objects_mask(image, background_level)

        max_value = np.max(image)
        if np.sum(working_mask) > 0:
            max_value = np.max(image[working_mask])

        # plt.imshow(working_mask, origin="lower")
        # plt.show()

        new_y, new_x = np.where(image == max_value)
        new_x = new_x[0] + closest_x - size
        new_y = new_y[0] + closest_y - size
        return new_x, new_y

    def reset_object_coords(self):
        """Recalculate the object coordinates."""
        for idx, _object in enumerate(self.obj_list):
            x, y = _object.xcoord, _object.ycoord
            size = self.max_radius
            image = self.image[y - size : y + size, x - size : x + size]
            # median = np.median(image)
            # std = np.std(image)
            # plt.imshow(
            #     image, origin="lower", vmax=median + 3 * std, vmin=median - 3 * std
            # )
            # plt.show()

            sky, sky_err = self._calc_background_level(image)
            background_level = sky + 4 * sky_err  # TODO: check this
            working_mask = self._create_objects_mask(image, background_level)
            closest_x, closest_y = x, y
            if np.sum(working_mask) > 0:
                closest_x, closest_y = self._find_closest_object(working_mask)
                closest_x += x - size
                closest_y += y - size

            new_x, new_y = self._find_coords_max_pixel(closest_x, closest_y)
            ra, dec = self._convert_pixels_to_coords(new_x, new_y)
            new_x += 1
            new_y += 1
            _object.xcoord, _object.ycoord = new_x, new_y
            _object.ra, _object.dec = ra, dec
            _object.sky_photons = sky
            self.obj_list[idx] = _object
        return

    def calc_psf_radius(self):
        """Calculate FWHM of the object

        Parameters
        ----------
        size : int, optional
            size of the box in which the object will be evaluated, by default 20.

        Returns
        -------
        float
            FWHM calculated for the object.
        """
        for idx, _object in enumerate(self.obj_list):
            x, y, sky_photons = _object.xcoord, _object.ycoord, _object.sky_photons
            r = self.max_radius
            img_data = self.image[y - r : y + r, x - r : x + r] - sky_photons
            # median = np.median(self.image)
            # std = np.std(self.image)
            # plt.imshow(
            #     img_data, origin="lower", vmax=median + 3 * std, vmin=median - 3 * std
            # )
            # plt.show()

            light_profile = np.take(img_data, r - 1, axis=0)
            half_max = np.max(light_profile) / 2
            n = len(light_profile)
            x = np.linspace(0, n - 1, n)
            spline = UnivariateSpline(x, light_profile - half_max, s=None)
            # plt.plot(spline(x))
            # plt.show()

            roots = spline.roots()
            idx_max_val = np.argmax(spline(x))
            tmp = np.abs(roots - idx_max_val)
            idx_r_1 = np.argmin(tmp)
            tmp[idx_r_1] = 1e10
            idx_r_2 = np.argmin(tmp)
            # plt.plot(spline(x))
            # plt.plot(roots[idx_r_2], 0, "bo")
            # plt.plot(roots[idx_r_1], 0, "bo")
            # plt.show()

            fwhm = np.abs(roots[idx_r_2] - roots[idx_r_1])
            self.obj_list[idx].psf_radius = 3 * fwhm

        return

    def _create_sky_mask(self, xcoord, ycoord, psf_radius):
        working_mask = np.ones(self.image_shape, bool)
        ym, xm = np.indices(self.image_shape, dtype="float32")
        r = np.sqrt((xm - xcoord) ** 2 + (ym - ycoord) ** 2)
        mask = (r > 2 * psf_radius) * (r < 3 * psf_radius) * working_mask
        return mask

    def calc_sky_photons(self):
        """Calculate the number of photons of the sky"""
        for idx, _object in enumerate(self.obj_list):
            mask = self._create_sky_mask(
                _object.xcoord, _object.ycoord, _object.psf_radius
            )
            sky = self.image[np.where(mask)]
            self.obj_list[idx].sky_photons = np.median(sky)

    def _create_psf_mask(self, xcoord, ycoord, psf_radius):
        working_mask = np.ones(self.image_shape, bool)
        ym, xm = np.indices(self.image_shape, dtype="float32")
        r = np.sqrt((xm - xcoord) ** 2 + (ym - ycoord) ** 2)
        mask = (r < psf_radius) * working_mask
        return mask

    def calc_psf_photons(self):
        """Calculate the number of photons of the object."""
        for idx, _object in enumerate(self.obj_list):
            mask = self._create_psf_mask(
                _object.xcoord, _object.ycoord, _object.psf_radius
            )
            star = self.image[np.where(mask)]
            star_photons = np.sum(star - _object.sky_photons)
            if star_photons < 0:
                star_photons = 0
            self.obj_list[idx].star_photons = star_photons
            self.obj_list[idx].star_err = np.sqrt(
                star_photons + _object.sky_photons * star.shape[0]
            )  # TODO: add the read noise error
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
    psf_radius: float = 0
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
