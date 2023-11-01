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
from sys import exit


class FITS_files_manager:
    """This class is a manager to deal with groups of FITS files."""

    def __init__(self, dir_path: str, file_name_tag: str = ".fits"):
        self.file_name_tag = file_name_tag
        if not os.path.isdir(dir_path):
            raise FileNotFoundError(dir_path)
        self.dir_path = dir_path
        self._create_FITS_objs()

        return

    def _create_FITS_objs(self):
        self.cam_files = {"cam3": [], "cam4": []}
        list_dir = [f for f in os.listdir(self.dir_path) if self.file_name_tag in f]
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

    def get_images_by_rotor_position(self, rpos: int) -> list:
        """Get set of images by the rotor positions.

        Parameters
        ----------
        rpos : int
            current rotor position.

        Returns
        -------
        list
            list of FITS_file objects.
        """
        current_rpos = {}
        for cam in [3, 4]:
            current_rpos[f"cam{cam}"] = [
                obj for obj in self.cam_files[f"cam{cam}"] if obj.rot_pos == rpos
            ]
        return current_rpos

    def shift_images(self, obj_coords_file: str):
        """Shift the images.

        Shift the images based on the coordinates of the object found in the 'obj_coords_file'.

        Parameters
        ----------
        obj_coords_file : str
            A csv file with the x and y coordiantes of the object over the image series.
        """
        shifts = self._get_shifts(obj_coords_file)
        dest_path = os.path.join(self.dir_path, "..", "shifted_images")
        if not os.path.exists(dest_path):
            os.mkdir(dest_path)

        for cam, ffiles in self.cam_files.items():
            for idx, ffile in enumerate(ffiles):
                x_shift = shifts[f"{cam}_x"][idx]
                y_shift = shifts[f"{cam}_y"][idx]
                file_name = os.path.join(self.dir_path, ffile.name)
                data, hdr = fits.getdata(file_name, header=True)
                data = self._shift_image(data, x_shift, y_shift)
                file_name = os.path.join(dest_path, ffile.name)
                fits.writeto(file_name, data, hdr, overwrite=True)
        self.dir_path = dest_path

        return

    def combine_images_by_run(self, dest_path: str):
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
            for _, ffiles in current_run.items():
                images = []
                for ffile in ffiles:
                    file_name = os.path.join(self.dir_path, ffile.name)
                    data, hdr = fits.getdata(file_name, header=True)
                    images.append(data)
                file_name = os.path.join(dest_path, ffile.name)
                median = np.median(images, axis=0)
                hdr["expnum"] = 0
                fits.writeto(file_name, median, hdr, overwrite=True)
        return

    def combine_images_by_rotor_position(self, dest_path: str, nruns=None):
        """Combine a set of images of the same rotor position.

        Parameters
        ----------
        dest_path : str
            destination path

        nruns : int, optional
            Number of runs to be combined. The default is None.

        use_moptp_name : bool, optional
            If True, use the MOPTOP name for the images. The default is False.
        """

        rotor_positions = self._get_rotor_positions()
        for rpos in rotor_positions:
            current_rpos = self.get_images_by_rotor_position(rpos)
            for ffiles in current_rpos.values():
                images = []
                for idx1, _tuple in enumerate(zip(*[iter(ffiles)] * nruns)):
                    for idx2, ffile in enumerate(_tuple):
                        idx = idx2 + idx1 * nruns
                        file_name = os.path.join(self.dir_path, ffile.name)
                        data, hdr = fits.getdata(file_name, header=True)
                        images.append(data)

                    file_name = os.path.join(dest_path, ffile.name)
                    median = np.mean(images, axis=0)
                    hdr["runnum"] = 0
                    fits.writeto(file_name, median, hdr, overwrite=True)

    @staticmethod
    def _get_shifts(shifts_file):
        df = pd.read_csv(shifts_file)
        shifts = df.drop(df.columns[[0, 1]], axis=1).to_dict(orient="list")

        for key, value in shifts.items():
            shifts[key] = np.asarray(shifts[key]) - value[0]

        return shifts

    @staticmethod
    def _shift_image(image, x_shift, y_shift):
        # TODO: implement the variable new_size
        xsize, ysize = image.shape
        x, y = xsize // 2 + x_shift, ysize // 2 + y_shift
        new_size = 480
        image = image[y - new_size : y + new_size + 1, x - new_size : x + new_size + 1]

        return image

    # @staticmethod
    # def _get_moptop_file_name(hdr):
    #     return hdr["EXPID"][:-1] + "1.fits"

    def _get_rotor_positions(self):
        rpos_numbers = [obj.rot_pos for obj in self.cam_files["cam3"]]
        counter = collections.Counter(rpos_numbers).items()
        imgs_per_rotor_position = [n for _, n in counter]
        if sum(imgs_per_rotor_position) % 16 != 0:
            raise ValueError(
                "There are not the same number of images per rotor position."
            )
        return [pos for pos, _ in counter]


@dataclass
class FITS_file:
    """This class keeps the header information needed to deal with the moptop polarimetry."""

    name: str
    mjd: float
    run_num: int
    rot_pos: int

    def __lt__(self, other):
        if isinstance(other, FITS_file):
            return self.mjd < other.mjd
