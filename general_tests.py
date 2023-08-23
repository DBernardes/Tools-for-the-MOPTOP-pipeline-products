from photometry import Photometry

#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


import os
import matplotlib.pyplot as plt
import numpy as np

import astropy.io.fits as fits

star_name = "GRB 230818A"
experiment = "all data"
src_path = os.path.join(
    "..",
    "Pol charact MOPTOP",
    "Scientific objects",
    star_name,
    experiment,
    star_name,
)

RA, DEC = "19:03:40", "+40:52:35"
size = 30
file = os.path.join(src_path, "3_e_20230818_6_3_4_1.fits")
phot = Photometry(file, RA, DEC)
print(phot.xcoord, phot.ycoord)
phot.recenter_object()
print(phot.xcoord, phot.ycoord)
phot.calc_psf_radius()
print(phot.psf_radius)
phot.calc_sky_photons()
print(phot.sky_photons)
star_photons = phot.calc_psf_photons()
print(star_photons)
