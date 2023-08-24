# This is the photometry class

import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import UnivariateSpline
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from dataclasses import dataclass
from pandas import DataFrame


class Photometry:
    """Photometry class"""

    def __init__(self, file: str, objects: DataFrame, max_radius: int = 30) -> None:
        """Initialize the class

        Parameters
        ----------
        file : str
            FITS file name.
        objects: DataFrame
            A pandas Dataframe of tuples with the name, right ascension, and declination of the objects.
        max_radius: int
            maximum radius, in pixels, in which the object can be found. Default to 30.
        """
        self.file = file
        self.max_radius = max_radius // 2
        self.image, self.header = fits.getdata(file, header=True)
        self.image_shape = self.image.shape
        self.obj_list = []
        for _object in objects.itertuples(name=None, index=False):
            name, ra, dec = _object
            xcoord, ycoord = self._convert_coords_to_pixel(ra, dec)
            self.obj_list.append(Object(name, xcoord, ycoord))
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

    def reset_object_coords(self):
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
                light_profile = np.take(img_data, r - 1, axis=0)
                max_star_flux = np.max(img_data)
                half_max = max_star_flux / 2
                n = len(light_profile)
                x = np.linspace(0, n, n)
                spline = UnivariateSpline(x, light_profile - half_max, s=0)
                r1, r2 = spline.roots()
                fwhm = r2 - r1
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
            _, xcoord, ycoord, psf_radius, *_ = _object.get_info()
            mask = self._create_sky_mask(xcoord, ycoord, psf_radius)
            self.obj_list[idx].sky_photons = np.median(self.image[np.where(mask)])

    def _create_psf_mask(self, xcoord, ycoord, psf_radius):
        working_mask = np.ones(self.image_shape, bool)
        ym, xm = np.indices(self.image_shape, dtype="float32")
        r = np.sqrt((xm - xcoord) ** 2 + (ym - ycoord) ** 2)
        mask = (r < psf_radius) * working_mask
        return mask

    def calc_psf_photons(self):
        """Calculate the number of photons of the object."""
        for idx, _object in enumerate(self.obj_list):
            _, xcoord, ycoord, psf_radius, *_ = _object.get_info()
            mask = self._create_psf_mask(xcoord, ycoord, psf_radius)
            self.obj_list[idx].star_photons = np.sum(
                self.image[np.where(mask)] - _object.sky_photons
            )
        return

    def get_mjd(self) -> float:
        """Get MJD of the image header

        Returns
        -------
        float
            Modified Julian Day (MJD)
        """
        return self.header["MJD"]


@dataclass
class Object:
    """Class to keep the photometru information related to an astronomic object."""

    name: str
    xcoord: str
    ycoord: str
    psf_radius: float = 0
    sky_photons: float = 0
    star_photons: float = 0

    def get_info(self):
        return (
            self.name,
            self.xcoord,
            self.ycoord,
            self.psf_radius,
            self.sky_photons,
            self.star_photons,
        )
