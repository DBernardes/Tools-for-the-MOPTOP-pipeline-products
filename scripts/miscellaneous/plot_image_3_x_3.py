import os
import matplotlib.pyplot as plt

xsize = 2048
ysize = 2048
div = 8
sep = 2
stepx = xsize // div
stepy = ysize // div
plate_scale = xsize / 420

for y in range(stepy, ysize, sep * stepy):
    for x in range(stepx, xsize, sep * stepx):
        plt.plot(x, y, "bo", alpha=0.5)
        strx = (x - stepx) / plate_scale
        stry = (y - stepy) / plate_scale
        plt.annotate(f"({strx:.2f},{stry:.2f})", (x, y + 30), fontsize=10, ha="center")


plt.xlim(0, xsize)
plt.ylim(0, ysize)
plt.show()
