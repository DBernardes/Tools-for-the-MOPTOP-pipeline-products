# This is the photometry class

import numpy as np
import astropy.io.fits as fits
from scipy.interpolate import UnivariateSpline
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u


class Photometry:
    """Photometry class"""

    def __init__(self, file: str, ra: int, dec: int, max_radius: int = 30) -> None:
        """Initialize the class

        Parameters
        ----------
        file : str
            FITS file name.
        ra : int
            right ascension of the object, in 'HH:MM:SS.ss'
        dec : int
            declination of the object, in 'HH:MM:SS.ss'
        max_radius:
            maximum radius, in pixels, in which the object can be found.
        """
        self.file = file
        self.max_radius = max_radius // 2
        self.image, self.header = fits.getdata(file, header=True)
        self.image_shape = self.image.shape
        self._convert_coords_to_pixel(ra, dec)
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
        self.xcoord = int(x) + 1
        self.ycoord = int(y) + 1
        return

    def recenter_object(self):
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
        x, y = self.xcoord, self.ycoord
        size = self.max_radius
        max_value = self._get_max_pixel_value()
        image = self.image[y - size : y + size, x - size : x + size]
        new_y, new_x = np.where(image == max_value)
        new_x += 1 + x - size
        new_y += 1 + y - size
        self.xcoord, self.ycoord = new_x[0], new_y[0]
        return

    def _get_max_pixel_value(self) -> float:
        """Get the value of the pixels inside a box

        Parameters
        ----------
        size : int, optional
            size of the box in which the max value will be seeked, by default 20.

        Returns
        -------
        float
            value of the pixels, according to the calculation method.
        """
        pixels_val = np.max(
            self.image[
                self.ycoord - self.max_radius : self.ycoord + self.max_radius,
                self.xcoord - self.max_radius : self.xcoord + self.max_radius,
            ]
        )
        return pixels_val

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
        x, y = self.xcoord, self.ycoord
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
        self.psf_radius = 3 * fwhm
        return

    def _create_sky_mask(self):
        working_mask = np.ones(self.image_shape, bool)
        ym, xm = np.indices(self.image_shape, dtype="float32")
        r = np.sqrt((xm - self.xcoord) ** 2 + (ym - self.ycoord) ** 2)
        mask = (r > 2 * self.psf_radius) * (r < 3 * self.psf_radius) * working_mask
        return mask

    def calc_sky_photons(self):
        """Calculate the number of photons of the sky"""
        mask = self._create_sky_mask()
        self.sky_photons = np.median(self.image[np.where(mask)])

    def _create_psf_mask(self):
        working_mask = np.ones(self.image_shape, bool)
        ym, xm = np.indices(self.image_shape, dtype="float32")
        r = np.sqrt((xm - self.xcoord) ** 2 + (ym - self.ycoord) ** 2)
        mask = (r < self.psf_radius) * working_mask
        return mask

    def calc_psf_photons(self) -> float:
        """Calculate the number of photons of the object.

        Returns
        -------
        float
            Number of photons of the object.
        """
        mask = self._create_psf_mask()
        star_photons = np.sum(self.image[np.where(mask)] - self.sky_photons)
        return star_photons
