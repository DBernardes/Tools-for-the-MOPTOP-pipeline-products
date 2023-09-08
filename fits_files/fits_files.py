#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"

import numpy as np
import astropy.io.fits as fits
from dataclasses import dataclass
import os
import collections
import pandas as pd
import matplotlib.pyplot as plt


class FITS_files_manager:
    """This class is a manager to deal with groups of FITS files."""

    def __init__(self, dir_path: str, file_name_tag: str = ".fits"):
        if not os.path.isdir(dir_path):
            raise FileNotFoundError(dir_path)
        self.dir_path = dir_path
        self._create_FITS_objs(file_name_tag)

        return

    def _create_FITS_objs(self, file_name_tag):
        self.cam_files = {"cam3": [], "cam4": []}
        list_dir = [f for f in os.listdir(self.dir_path) if file_name_tag in f]
        for file in list_dir:
            file_path = os.path.join(self.dir_path, file)
            hdr = fits.getheader(file_path)
            ffile = FITS_file(file, hdr["MJD"], hdr["RUNNUM"], hdr["EXPNUM"])

            self.cam_files[f"cam{hdr['EXPID'][0]}"].append(ffile)
        self.cam_files["cam3"].sort()
        self.cam_files["cam4"].sort()

    def get_images_by_run(self, run: int) -> list:
        """Get set of images by the run number.

        Parameters
        ----------
        run : int
            current run number.

        Returns
        -------
        list
            list of FITS_file objects.
        """
        current_run = {}
        for cam in [3, 4]:
            current_run[f"cam{cam}"] = [
                obj for obj in self.cam_files[f"cam{cam}"] if obj.run_num == run
            ]
        return current_run

    def combine_images_by_run(self, dest_path: str, shifts_file: str = ""):
        """Combine a set of images of the same run.

        Parameters
        ----------
        dest_path : str
            destination path;
        """
        run_numbers = [obj.run_num for obj in self.cam_files["cam3"]]
        run_numbers = [item for item, _ in collections.Counter(run_numbers).items()]
        for run in run_numbers:
            current_run = self.get_images_by_run(run)
            for cam, ffiles in current_run.items():
                images = []
                shifts = self._get_shifts(run, shifts_file)
                for idx, ffile in enumerate(ffiles):
                    file_name = os.path.join(self.dir_path, ffile.name)
                    data, hdr = fits.getdata(file_name, header=True)
                    if shifts_file != "":
                        x_shift, y_shift = (
                            shifts[f"{cam}_x"][idx],
                            shifts[f"{cam}_y"][idx],
                        )
                        data = self._shift_image(data, x_shift, y_shift)
                    images.append(data)
                file_name = os.path.join(dest_path, f"{cam[-1]}_e_run{run}.fits")
                median = np.median(images, axis=0)
                hdr["expnum"] = 0
                fits.writeto(file_name, median, hdr, overwrite=True)
        return

    @staticmethod
    def _get_shifts(run, shifts_file):
        if shifts_file == "":
            return []
        else:
            df = pd.read_csv(shifts_file)
            rows = df.loc[df["run_num"] == run]
            shifts = {}
            for name, *val in rows.transpose().itertuples(name=None):
                shifts[name] = np.asarray(val)

            shifts["cam3_x"] -= shifts["cam3_x"][0]
            shifts["cam3_y"] -= shifts["cam3_y"][0]
            shifts["cam4_x"] -= shifts["cam4_x"][0]
            shifts["cam4_y"] -= shifts["cam4_y"][0]
            return shifts

    @staticmethod
    def _shift_image(image, x_shift, y_shift):
        xsize, ysize = image.shape
        x, y = xsize // 2 + x_shift, ysize // 2 + y_shift

        image = image[y - 500 : y + 500 + 1, x - 500 : x + 500 + 1]

        return image


@dataclass
class FITS_file:
    """This class keeps the header information needed to deal with the moptop polarimetry."""

    name: str
    mjd: float
    run_num: int
    exp_num: int

    def __lt__(self, other):
        if isinstance(other, FITS_file):
            return self.mjd < other.mjd
