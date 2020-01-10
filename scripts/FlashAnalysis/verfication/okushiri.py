#!/usr/bin/python3

# Verification of the okushiri scenario in [Vater, 2019]

import sys, core
from signal import *
from pysamoaxdmf.reader import Reader, Parameter
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy import interpolate
from os.path import basename, dirname
import pandas as pd
from multiprocessing import Pool
from time import time
import itertools
import interp

gauges = [5, 7, 9]
ts = [15.0, 15.5, 16.0, 16.5, 17.0]

levelsh = np.linspace(0.0, 0.04, 100)
cmapb = matplotlib.colors.LinearSegmentedColormap.from_list("", ["limegreen","peru","white"])

def sfplot(argv, name, filename):
    # Use getopt to parse argv
    inputfile, _, outputdir, _, _ = core.args(argv)
    print('Reading XMF file: ' + inputfile)
    print('Writing output surface graphs to: ' + outputdir)

    # Read data file and load steps
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    print(str(xdmf) + "\n")
    sxs = list(map(xdmf.step_near_time, ts))
    for sx, t in zip(sxs, ts):  
        xdmf.steps[sx].load()
        print(str(xdmf.steps[sx]))
        print('   index: ' + str(sx) + ', time deviation: ' + str(xdmf.steps[sx].time - t) + 's')

    # Convert data into plottable format
    pts = list(map(lambda sx: xdmf.steps[sx].cells_x.reshape(-1, xdmf.steps[sx].cells_x.shape[-1]), sxs))
    bs = list(map(lambda sx: xdmf.steps[sx].cells_b.reshape(-1), sxs))
    hs = list(map(lambda sx: xdmf.steps[sx].cells_h.reshape(-1), sxs))
    hbs = list(map(lambda pack: [a + b - 0.32 if b < 0.32 else a for a, b in zip(pack[0], pack[1])], zip(hs, bs)))

    # Plot results
    plt.rcParams['figure.figsize'] = [20, 15]
    fig, subplots = plt.subplots(2, 3)
    fig.delaxes(subplots[-1, -1])
    fig.suptitle(name + "\n" + 'Surface plot')
    for i, subplot in enumerate(subplots.flat):
        if i < len(ts):
            print('Plot at t = ' + str(ts[i]))
            xsf, zsf = zip(*(pts[i]))
            subplot.set_title('t = ' + str(ts[i]))
            subplot.set(xlim=(4.6, 5.448), ylim=(1.4, 2.4))
            subplot.tricontourf(xsf, zsf, bs[i], cmap=cmapb, levels=10)
            subplot.tricontour(xsf, zsf, bs[i], levels=10, linewidths=0.5, colors='k')
            triag = tri.Triangulation(xsf, zsf)
            triag.set_mask([np.any(hs[i][tri] < 0.0001) for tri in triag.triangles])
            tpc = subplot.tricontourf(triag, hs[i], cmap="coolwarm", levels=levelsh, extend='max')
            for c in tpc.collections:
                c.set_edgecolor("face")
    fig.colorbar(tpc, ax=subplots.ravel().tolist(), shrink=0.5)
    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(20, 15)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')

if __name__ == '__main__':
    signal(SIGINT, core.siginth)
    _, _, _, scenario, _ = core.args(sys.argv[1:])
    if (scenario == 'surface'):
        sfplot(sys.argv[1:], 'Runup onto a complex three-dimensional beach', 'surface')
    else:
        interp.csvcompexp(sys.argv[1:], 'Runup onto a complex three-dimensional beach', 'csvcomp', gauges, \
            'gauges.csv', 0.01, 0.0, 0.0, 0.0, 40.0, -0.02, 0.05)
