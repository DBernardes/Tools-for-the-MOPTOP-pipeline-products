from tools import track_obj_over_images
import os
import matplotlib.pyplot as plt

xsize = 1024
ysize = 1024
div = 8
sep = 3
stepx = xsize//div
stepy = ysize//div

for y in range(stepy, ysize, sep*stepy):
    for x in range(stepx, xsize, sep*stepx):
        plt.plot(x, y, 'bo', alpha=0.5)
        plt.annotate(f'({x},{y})', (x,y+20), fontsize=10, ha='center')


plt.xlim(0, xsize)
plt.ylim(0, ysize)
plt.show()