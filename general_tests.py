import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

star_name = "GRB 230818A"
experiment = "all data"
src_path = os.path.join(
    "..", "Pol charact MOPTOP", "Scientific objects", star_name, experiment
)
csv_path = os.path.join(src_path, "csv", "second set")


csv_file = os.path.join(csv_path, "original.csv")
df = pd.read_csv(csv_file)
original = df["star_photons"]
csv_file = os.path.join(csv_path, "comparison_2.csv")
df = pd.read_csv(csv_file)
comp_1 = df["star_photons"]
csv_file = os.path.join(csv_path, "comparison_2.csv")
df = pd.read_csv(csv_file)
comp_2 = df["star_photons"]
csv_file = os.path.join(csv_path, "candidate_1.csv")
df = pd.read_csv(csv_file)
cand_1 = df["star_photons"]
csv_file = os.path.join(csv_path, "candidate_2.csv")
df = pd.read_csv(csv_file)
cand_2 = df["star_photons"]
csv_file = os.path.join(csv_path, "candidate_3.csv")
df = pd.read_csv(csv_file)
cand_3 = df["star_photons"]

# comp = comp_2 + comp_1
# cand_1 /= comp
# # cand_1 /= np.median(cand_1)
# cand_2 /= comp
# # cand_2 /= np.median(cand_2)
# cand_3 /= comp
# # cand_3 /= np.median(cand_3)

plt.plot(df["mjd"], comp_1, "bo", alpha=0.5)
plt.show()
