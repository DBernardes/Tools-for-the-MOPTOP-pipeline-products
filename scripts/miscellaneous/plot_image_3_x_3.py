import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

xsize = 1024
ysize = xsize
div = 8
sep = 2
stepx = xsize // div
stepy = ysize // div
plate_scale = xsize / 420
initial_stepx = (0.5 * xsize + stepx) / plate_scale
initial_stepy = (0.5 * ysize - stepy) / plate_scale

fig, ax = plt.subplots()

for y in range(-stepy, ysize + 2 * stepy, sep * stepy):
    for x in range(-stepx, xsize + 2 * stepx, sep * stepx):
        ax.plot(x, y, "bo", alpha=0.5)
        strx = x + stepx
        stry = y + stepy
        ax.annotate(f"({strx:.2f},{stry:.2f})", (x, y + 30), fontsize=10, ha="center")

# Create a Rectangle patch
rect = patches.Rectangle(
    (0, 0),
    xsize,
    ysize,
    linewidth=1,
    edgecolor="r",
    facecolor="none",
)

plt.title(f"Initial step: -{initial_stepx:.2f},-{initial_stepy:.2f}")
ax.add_patch(rect)
plt.xlim(0, xsize)
plt.ylim(0, ysize)
# plt.title("Position of the object over the CCD frame")
plt.xlabel("X axis (pixels)")
plt.ylabel("Y axis (pixels)")
plt.show()
