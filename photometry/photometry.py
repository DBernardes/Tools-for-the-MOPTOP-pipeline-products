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
from math import sqrt


class Photometry:
    """Photometry class"""

    def __init__(self, file: str, objects: DataFrame, max_size: int = 10) -> None:
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
    ) -> None:
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

    def _calc_background_level(self, image: np.ndarray, nsigma: int = 5) -> float:
        median = np.median(image)
        std = np.median(np.abs(image - median))

        indexes = np.where(
            (median - nsigma * std < image) & (image < median + nsigma * std)
        )
        background_pixels = image[indexes]
        median = np.median(background_pixels)
        std = np.std(background_pixels)
        self.background_level = median + nsigma * std
        # print(self.background_level)

    def _create_objects_maks(self, image):
        working_mask = np.zeros(image.shape, bool)
        y_coords, x_coords = np.where(image > self.background_level)
        working_mask[(y_coords, x_coords)] = 1

        for x, y in zip(x_coords, y_coords):
            sub_image = working_mask[y - 1 : y + 2, x - 1 : x + 2]
            if np.sum(sub_image) == 1:
                working_mask[y, x] = 0
        self.working_mask = working_mask

        return

    def _find_closest_object(self):
        center_coord = self.working_mask.shape[0] // 2
        y_coords, x_coords = np.where(self.working_mask == 1)
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

        self._calc_background_level(image)
        working_mask = np.zeros(image.shape, bool)
        y_coords, x_coords = np.where(image > self.background_level)
        working_mask[(y_coords, x_coords)] = 1
        for x, y in zip(x_coords, y_coords):
            sub_image = working_mask[y - 1 : y + 2, x - 1 : x + 2]
            if np.sum(sub_image) == 1:
                working_mask[y, x] = 0

        plt.imshow(working_mask)
        plt.show()

        if np.sum(working_mask) == 0:
            max_value = np.max(image)
        else:
            max_value = np.max(image[working_mask])
        new_y, new_x = np.where(image == max_value)
        new_x = new_x[0] + closest_x - size
        new_y = new_y[0] + closest_y - size
        return new_x, new_y

    def reset_object_coords(self):
        """Recalculate the object coordinates."""
        for idx, _object in enumerate(self.obj_list):
            _, x, y, *_ = _object.get_info()
            size = self.max_radius
            image = self.image[y - size : y + size, x - size : x + size]
            self._calc_background_level(image)
            self._create_objects_maks(image)
            if np.sum(self.working_mask) == 0:
                closest_x, closest_y = x, y
            else:
                closest_x, closest_y = self._find_closest_object()
                closest_x += x - size
                closest_y += y - size
                # print(closest_x, closest_y, "\n")
            new_x, new_y = self._find_coords_max_pixel(closest_x, closest_y)
            new_x += 1
            new_y += 1
            self.obj_list[idx].xcoord, self.obj_list[idx].ycoord = new_x, new_y
        return

    def reset_object_coords_1(self):
        """Recalculate the object coordinates.

        Parameters
        ----------
        size : int, optional
            size of the box in which the object will be evaluated, by default 20.

        Returns
        -------
        tuple[int, int]
            x and y coordinates of the pixel.
        """
        for idx, _object in enumerate(self.obj_list):
            _, x, y, *_ = _object.get_info()
            size = self.max_radius
            image = self.image[y - size : y + size, x - size : x + size]
            max_value = np.max(image)
            new_y, new_x = np.where(image == max_value)
            new_x += 1 + x - size
            new_y += 1 + y - size
            self.obj_list[idx].xcoord, self.obj_list[idx].ycoord = new_x[0], new_y[0]
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
            try:
                _, x, y, *_ = _object.get_info()
                r = self.max_radius
                img_data = self.image[y - r : y + r, x - r : x + r]
                # plt.imshow(img_data)
                # plt.show()
                light_profile = np.take(img_data, r - 1, axis=0)
                max_star_flux = np.max(img_data)
                half_max = max_star_flux / 2
                n = len(light_profile)
                x = np.linspace(0, n, n)
                spline = UnivariateSpline(x, light_profile - half_max, s=0)
                roots = spline.roots()

                idx_max_val = np.argmax(spline(x))
                tmp = np.abs(roots - idx_max_val)
                idx_r_1 = np.argmin(tmp)
                tmp[idx_r_1] = 1e10
                idx_r_2 = np.argmin(tmp)
                # plt.plot(spline(x))
                # print(roots[idx_r_2], roots[idx_r_1])
                # plt.show()

                fwhm = np.abs(roots[idx_r_2] - roots[idx_r_1])
                self.obj_list[idx].psf_radius = 3 * fwhm
            except Exception:
                continue
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
            _, xcoord, ycoord, _, psf_radius, *_ = _object.get_info()
            mask = self._create_sky_mask(xcoord, ycoord, psf_radius)
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
            (
                _,
                xcoord,
                ycoord,
                _,
                psf_radius,
                sky_photons,
                *_,
            ) = _object.get_info()
            mask = self._create_psf_mask(xcoord, ycoord, psf_radius)
            star = self.image[np.where(mask)]
            star_photons = np.sum(star - sky_photons)
            self.obj_list[idx].star_photons = star_photons
            self.obj_list[idx].star_err = np.sqrt(
                star_photons + sky_photons * star.shape[0]
            )
        return


@dataclass
class Object:
    """Class to keep the photometry information related to an astronomic object."""

    name: str
    xcoord: str
    ycoord: str
    mjd: float
    psf_radius: float = 0
    sky_photons: float = 0
    star_photons: float = 0
    star_err: float = 0

    def get_info(self):
        return (
            self.name,
            self.xcoord,
            self.ycoord,
            self.mjd,
            self.psf_radius,
            self.sky_photons,
            self.star_photons,
            self.star_err,
        )
