# import os
# import pandas as pd
# from photometry import Photometry
# from scripts.tools import sort_files
# import astropy.io.fits as fits
# import matplotlib.pyplot as plt
# import numpy as np
# from photometry import Photometry

# star_name = "GRB 230818A"
# experiment = "all data"
# _set = "first"
# camera = 3
# src_path = os.path.join(
#     "..",
#     "Pol charact MOPTOP",
#     "Scientific objects",
#     star_name,
#     experiment,
#     star_name,
#     f"{_set} set",
# )


# csv_file = os.path.join(src_path, "..", f"objects coordinates.csv")
# df = pd.read_csv(csv_file)
# objects = {
#     "name": df["name"],
#     "ra": df[f"ra_{_set}_set_cam{camera}"],
#     "dec": df[f"dec_{_set}_set_cam{camera}"],
# }
# objects = pd.DataFrame.from_dict(objects)

# objects_photometry = {}
# for obj_name in objects["name"]:
#     objects_photometry[obj_name] = {
#         "mjd": [],
#         "xcoord": [],
#         "ycoord": [],
#         "star_photons": [],
#     }

# file = "3_e_20230818_5_2_1_1.fits"
# file_path = os.path.join(src_path, file)

# phot = Photometry(file_path, objects, 20)
# phot.reset_object_coords()
# phot.calc_psf_radius()
# phot.calc_sky_photons()
# phot.calc_psf_photons()
# mjd = phot.get_mjd()

# for _object in phot.obj_list:
#     print(repr(_object))


import scipy
import numpy as np
import matplotlib.pyplot as plt


def monoExp(x, m, t, b):
    return m * np.exp(-t * x) + b


xs = np.asarray(
    [
        60174.979717,
        60174.980688,
        60174.981614,
        60174.98254,
        60174.983466,
        60174.984392,
        60174.985318,
        60174.986244,
        60174.988095,
        60174.989021,
        60174.991799,
        60174.992725,
        60174.993651,
        60174.996429,
        60174.997355,
    ]
)
xs -= xs[0]
ys = np.asarray(
    [
        2.00386286,
        1.22964511,
        0.92545008,
        1.23599087,
        0.67159911,
        1.0,
        0.17513408,
        0.27170592,
        0.1546637,
        -0.00405001,
        0.10111231,
        0.35666107,
        0.19001122,
        0.13064874,
        0.1316298,
    ]
)
params, cv = scipy.optimize.curve_fit(monoExp, xs, ys)
m, t, b = params

print(params)
# plot the results
plt.plot(xs, ys, ".", label="data")
plt.plot(xs, monoExp(xs, m, t, b), "--", label="fitted")
plt.title("Fitted Exponential Curve")
plt.show()
