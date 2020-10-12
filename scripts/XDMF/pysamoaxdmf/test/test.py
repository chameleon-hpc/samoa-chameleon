#!/usr/bin/python3

import matplotlib.pyplot as plt
from pysamoaxdmf.reader import Reader
import numpy as np
from PIL import Image

from pysamoaxdmf import sampler

# inputfile = "/mnt/scratch/samoa/ader-dg/oscillating_lake/fat/unlimited/swe_20200928_140446.295_xdmf.xmf"
# inputfile = "/mnt/scratch/samoa/ader-dg/okushiri/unlimited/swe_20201005_111510.268_xdmf.xmf"
inputfile = "/home/christoph/samoa/samoa-ader-dg/output/okushiri/swe_20201012_033636.820_xdmf.xmf"


xdmf = Reader(inputfile)
num_steps = len(xdmf.steps)
print(xdmf)
sid = 20 # xdmf.step_near_time(0.01)
xdmf.steps[sid].load()
print(xdmf.steps[sid])
bnd = xdmf.steps[sid].bounds()
buffer = xdmf.steps[sid].sample_full(500.0)
xdmf.steps[sid].unload()

fig, ax = plt.subplots(figsize=(12, 12))
im = ax.imshow(buffer[:,:,3], cmap="coolwarm", aspect="equal", interpolation="bilinear", 
    origin="lower", extent=(bnd[0], bnd[2], bnd[1], bnd[3]))
fig.colorbar(im, ax=ax, shrink=0.5)

plt.savefig('test.png')
