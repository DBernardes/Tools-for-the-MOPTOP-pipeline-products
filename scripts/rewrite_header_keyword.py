from tools import rewrite_header_keyword, sort_files
import os
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u


star_name = "HD14069"
experiment = "several positions in image/20231117"
base_path = os.path.join(
    "..",
    "..",
    "Pol charact MOPTOP",
    "Low polarized stars",
    star_name,
    experiment,
    star_name,
)

# files = sort_files(base_path)
# file = files[0]
# file_before = os.path.join(base_path, file)
# hdr = fits.getheader(file_before)
# # ra = "2:16:27.548"
# # dec = "+07:36:48.13"
# for file in files:
#     file_after = os.path.join(base_path, file)
#     data = fits.getdata(file_after)
#     fits.writeto(file_after, data, hdr, overwrite=True)

# fits.setval(file_after, "CAT-RA", value=ra)
# fits.setval(file_after, "CAT-DEC", value=dec)  # .to_string(unit=u.deg, sep=":"))
# fits.setval(file_after, "RA", value=ra)
# fits.setval(file_after, "DEC", value=dec)  # .to_string(unit=u.deg, sep=":"))
# fits.setval(file_after, "WRA", value=ra)
# fits.setval(file_after, "WDEC", value=dec)  # .to_string(unit=u.deg, sep=":"))


i = -1
files = sort_files(base_path)
for idx_row in range(6):
    for idx_col in range(5):
        idx = idx_col + 6 * idx_row
        file = files[16 * idx]
        file_before = os.path.join(base_path, file)
        hdr = fits.getheader(file_before)
        ra, dec = hdr["CAT-RA"], hdr["CAT-DEC"]
        dec = Angle(dec, unit=u.deg) + i * Angle(f"105s")

        for file in files[(idx + 1) * 16 : (idx + 2) * 16]:
            file_after = os.path.join(base_path, file)
            fits.setval(file_after, "CAT-RA", value=ra)
            fits.setval(file_after, "CAT-DEC", value=dec.to_string(unit=u.deg, sep=":"))

    i *= -1
    hdr = fits.getheader(file_after)
    ra, dec = hdr["CAT-RA"], hdr["CAT-DEC"]
    ra = Angle(ra, unit=u.hour) - Angle(f"105s", unit=u.hour)
    for file in files[(idx + 5) * 16 : (idx + 6) * 16]:
        file_after = os.path.join(base_path, file)
        fits.setval(file_after, "CAT-DEC", value=dec)
        fits.setval(file_after, "CAT-RA", value=ra.to_string(unit=u.hour, sep=":"))
