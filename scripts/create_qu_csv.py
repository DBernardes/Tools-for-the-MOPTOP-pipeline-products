#!/usr/bin/env python

__author__ = "Denis Bernardes"
__copyright__ = "Copyright 2023, Liverpool John Moores University"


from tools import sort_qu_per_filter, low_polarized_stars, sigma_clipping, LIMIT_MJD
import os
import numpy as np
import pandas as pd


mean_qu_values = {"star": [], "filter": [], "q": [], "std_q": [], "u": [], "std_u": []}

for star in low_polarized_stars.keys():
    base_path = os.path.join(
        "..",
        "..",
        "Pol charact MOPTOP",
        "Low polarized stars",
        star,
        "all data",
        "reduced",
        star,
    )
    csv_file_name = os.path.join(base_path, "manipulated_data.csv")
    df = pd.read_csv(csv_file_name)
    df = df.loc[df["mjd"] > LIMIT_MJD]

    for filter in ["B", "V", "R", "I", "L"]:
        rows = df.loc[df["wave"] == f"MOP-{filter}"]
        if (
            filter == "B" and star != "BD+32 3739"
        ):  # this because we do not have enough data for the filter B for the other stars
            continue
        if rows.shape[0] == 0:
            continue
        q = rows["q_avg"]
        u = rows["u_avg"]
        q, u, *_ = sigma_clipping(q, u)

        meanq, meanu = np.mean(q), np.mean(u)
        stdq, stdu = np.std(q), np.std(u)
        mean_qu_values["star"].append(star)
        mean_qu_values["filter"].append(filter)
        mean_qu_values["q"].append(meanq)
        mean_qu_values["std_q"].append(stdq)
        mean_qu_values["u"].append(meanu)
        mean_qu_values["std_u"].append(stdu)

tmp_df = pd.DataFrame.from_dict(mean_qu_values)
tmp_dict = {"filter": [], "q": [], "std_q": [], "u": [], "std_u": []}
for filter in ["B", "V", "R", "I", "L"]:
    rows = tmp_df.loc[tmp_df["filter"] == filter]
    tmp_dict["filter"].append(filter)
    tmp_dict["q"].append(np.mean(rows["q"]))
    tmp_dict["std_q"].append(np.mean(rows["std_q"]))
    tmp_dict["u"].append(np.mean(rows["u"]))
    tmp_dict["std_u"].append(np.mean(rows["std_u"]))

csv_file = os.path.join(base_path, *4 * [".."], "mean_qu_values.csv")
df = pd.DataFrame.from_dict(tmp_dict)
pd.DataFrame.to_csv(df, csv_file, index=False)
