# # -*- coding: utf-8 -*-

import pytest
import os
from fits_files import FITS_files_manager, FITS_file
import astropy.io.fits as fits
import pandas as pd
import numpy as np
import collections


@pytest.fixture
def ffile():
    return FITS_file("test.fits", 0, 1, 1)


@pytest.fixture
def ffile_2():
    return FITS_file("test2.fits", 1, 2, 2)


def test_ffile_init(ffile):
    assert ffile.name == "test.fits"
    assert ffile.mjd == 0
    assert ffile.run_num == 1
    assert ffile.rot_pos == 1


def test_larger_than(ffile, ffile_2):
    assert ffile < ffile_2


# ----------------------------------------------------------------

dir_path = "./tests/images"


@pytest.fixture
def ffile_man():
    return FITS_files_manager(dir_path)


def test_ffile_man_init(ffile_man):
    assert ffile_man.dir_path == dir_path
    assert ffile_man.file_name_tag == ".fits"


def test_create_FITS_objects(ffile_man):
    cam_files = {"cam3": [], "cam4": []}
    ffile_man._create_FITS_objs()
    for file in os.listdir(dir_path):
        file_name = os.path.join(dir_path, file)
        hdr = fits.getheader(file_name)
        ffile = FITS_file(file, hdr["MJD"], hdr["RUNNUM"], hdr["EXPNUM"])
        cam_files[f"cam{hdr['EXPID'][0]}"].append(ffile)
    cam_files["cam3"].sort()
    cam_files["cam4"].sort()

    assert ffile_man.cam_files == cam_files


def test_get_images_by_run(ffile_man):
    run_num = 2
    cam_files = {"cam3": [], "cam4": []}
    current_run = ffile_man.get_images_by_run(run_num)
    for file in os.listdir(dir_path):
        file_name = os.path.join(dir_path, file)
        hdr = fits.getheader(file_name)
        ffile = FITS_file(file, hdr["MJD"], hdr["RUNNUM"], hdr["EXPNUM"])
        if hdr["RUNNUM"] == run_num:
            cam_files[f"cam{hdr['EXPID'][0]}"].append(ffile)
    cam_files["cam3"].sort()
    cam_files["cam4"].sort()
    assert current_run == cam_files


def test_get_images_by_rotor_pos(ffile_man):
    rpos = 2
    cam_files = {"cam3": [], "cam4": []}
    current_run = ffile_man.get_images_by_rotor_position(rpos)
    for file in os.listdir(dir_path):
        file_name = os.path.join(dir_path, file)
        hdr = fits.getheader(file_name)
        ffile = FITS_file(file, hdr["MJD"], hdr["RUNNUM"], hdr["EXPNUM"])
        if hdr["EXPNUM"] == rpos:
            cam_files[f"cam{hdr['EXPID'][0]}"].append(ffile)
    cam_files["cam3"].sort()
    cam_files["cam4"].sort()
    assert current_run == cam_files


obj_coords_file = "./tests/setup/coords_over_series.csv"


def test_get_shifts(ffile_man):
    shifts = ffile_man._get_shifts(obj_coords_file)
    df = pd.read_csv(obj_coords_file)
    new_shifts = df.drop(df.columns[[0, 1]], axis=1).to_dict(orient="list")
    for key, value in new_shifts.items():
        shifts[key].all() == (np.asarray(new_shifts[key]) - value[0]).all()


def test_shift_unique_image(ffile_man):
    image_path = os.path.join(dir_path, os.listdir(dir_path)[0])
    image = fits.getdata(image_path)
    shifts = ffile_man._get_shifts(obj_coords_file)

    x_shift, y_shift = shifts["cam3_x"][0], shifts["cam3_y"][0]
    ffile_man._shift_image(image, x_shift, y_shift)

    xsize, ysize = image.shape
    x, y = xsize // 2 + x_shift, ysize // 2 + y_shift
    new_size = 480
    image = image[y - new_size : y + new_size + 1, x - new_size : x + new_size + 1]

    return


def test_shift_images(ffile_man):
    ffile_man.shift_images(obj_coords_file)

    shifts = ffile_man._get_shifts(obj_coords_file)
    for cam, ffiles in ffile_man.cam_files.items():
        for idx, ffile in enumerate(ffiles):
            x_shift = shifts[f"{cam}_x"][idx]
            y_shift = shifts[f"{cam}_y"][idx]
            file_name = os.path.join(dir_path, ffile.name)
            data = fits.getdata(file_name)
            data = ffile_man._shift_image(data, x_shift, y_shift)

            fits_file_name = os.path.join(dir_path, "..", "shifted_images", ffile.name)
            fits_file = fits.getdata(fits_file_name)
            assert fits_file.all() == data.all()


def test_combine_images_by_run(ffile_man):
    new_path = os.path.join("tests", "shifted_images")
    dest_path = os.path.join(new_path, "..", "tmp")
    ffile_man.combine_images_by_run(dest_path)

    run_numbers = [1, 2]
    for run in run_numbers:
        current_run = ffile_man.get_images_by_run(run)
        for ffiles in current_run.values():
            images = []
            for ffile in ffiles:
                file_name = os.path.join(new_path, ffile.name)
                data = fits.getdata(file_name)
                images.append(data)
            median = np.median(images, axis=0)

            fits_file_name = os.path.join(dest_path, ffile.name)
            fits_file = fits.getdata(fits_file_name)
            assert fits_file.all() == median.all()
            os.remove(fits_file_name)


def test_verify_rotor_positions(ffile_man):
    ffile_man._verify_rotor_positions()


def test_combine_images_by_rotor_positions(ffile_man):
    nruns = 2
    dest_path = os.path.join(dir_path, "..", "tmp")
    ffile_man.combine_images_by_rotor_position(dest_path, nruns)

    for rpos in np.linspace(1, 16, 16, dtype=int):
        current_rpos = ffile_man.get_images_by_rotor_position(rpos)
        for ffiles in current_rpos.values():
            for _tuple in zip(*[iter(ffiles)] * nruns):
                images = []
                for ffile in _tuple:
                    file_name = os.path.join(dir_path, ffile.name)
                    data = fits.getdata(file_name)
                    images.append(data)
                median = np.mean(images, axis=0)

                fits_file_name = os.path.join(dest_path, ffile.name)
                fits_file = fits.getdata(fits_file_name)
                assert fits_file.all() == median.all()
                os.remove(fits_file_name)
