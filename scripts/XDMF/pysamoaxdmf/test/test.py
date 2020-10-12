#!/usr/bin/python3

import matplotlib.pyplot as plt
from pysamoaxdmf.reader import Reader
import numpy as np
from PIL import Image

# inputfile = "/mnt/scratch/samoa/ader-dg/oscillating_lake/fat/unlimited/swe_20200928_140446.295_xdmf.xmf"

# inputfile = "/mnt/scratch/samoa/ader-dg/okushiri/unlimited/swe_20201005_111510.268_xdmf.xmf"
inputfile = "/home/christoph/samoa/samoa-ader-dg/output/okushiri/swe_20201012_033636.820_xdmf.xmf"

xdmf = Reader(inputfile)
num_steps = len(xdmf.steps)
print(xdmf)
sid = 20 # xdmf.step_near_time(0.01)
xdmf.steps[sid].load()
print(xdmf.steps[sid])
# print(xdmf.steps[sid].sample((5.3, 2.3))) # 5.2

bnd = xdmf.steps[sid].bounds()
width = int((bnd[2] - bnd[0]) * 500.0)
height = int((bnd[3] - bnd[1]) * 500.0)
print(str(width) + "x" + str(height))
buffer = np.zeros((height, width))
for y in range(0, height):
    ny = bnd[1] + ((float(y) / height) * (bnd[3] - bnd[1]))
    print(y)
    for x in range(0, width):
        nx = bnd[0] + ((float(x) / width) * (bnd[2] - bnd[0]))
        s = xdmf.steps[sid].sample((nx, ny))
        buffer[y, x] = s[2]


xdmf.steps[sid].unload()

fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(buffer, cmap="coolwarm", aspect="equal", interpolation="bilinear", 
    origin="lower", extent=(bnd[0], bnd[2], bnd[1], bnd[3]))
fig.colorbar(im, ax=ax, shrink=0.5)

plt.savefig('test.png')

