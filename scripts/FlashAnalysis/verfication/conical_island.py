#!/usr/bin/python3

# Verification of the conical island scenario in [Vater, 2019]

import sys, core
from signal import *
from pysamoaxdmf.reader import Reader, Parameter
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from os.path import basename, dirname
import pandas as pd
from multiprocessing import Pool
from time import time
import itertools
import interp

gauges = [6, 9, 16, 22]

def sfplot(argv, name, filename, ts, zmin, zmax):
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
    plt.rcParams['figure.figsize'] = [16, 8]
    fig, subplots = plt.subplots(1, 2)
    fig.suptitle(name + "\n" + 'Surface plot')
    for i, subplot in enumerate(subplots):
        print('Plot at t = ' + str(ts[i]))
        xsf, zsf = zip(*(pts[i]))
        subplot.set_title('t = ' + str(ts[i]))
        subplot.set_aspect('equal')
        levels = np.linspace(zmin, zmax, 50)
        tpc = subplot.tricontourf(xsf, zsf, hbs[i], cmap="jet", levels=levels, extend='both')
        for c in tpc.collections:
            c.set_edgecolor("face")
        fig.colorbar(tpc, ax=subplot, shrink=0.5)
        c1 = plt.Circle((12.96, 13.8), 2.2 / 2.0, color='black', fill=False, lw=1)
        c2 = plt.Circle((12.96, 13.8), 4.64 / 2.0, color='black', fill=False, lw=1)
        c3 = plt.Circle((12.96, 13.8), 7.2 / 2.0, color='black', fill=False, lw=1)
        subplot.add_artist(c1)
        subplot.add_artist(c2)
        subplot.add_artist(c3)
    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(16, 8)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')


if __name__ == '__main__':
    signal(SIGINT, core.siginth)
    _, _, _, scenario, _ = core.args(sys.argv[1:])
    scen = scenario.split('_')
    if(scen[0] == 'surface'):
        ts = [13, 16]
        zmin = -0.003
        zmax = 0.024
        if(scen[1] == 'c'): 
            ts = [11, 14]
            zmin = -0.02
            zmax = 0.06
        sfplot(sys.argv[1:], 'Flow around a conical island, case ' + scen[1], 'surface', ts, zmin, zmax)
    else:
        ymin = -0.015
        ymax = 0.025
        if(scen[1] == 'c'): 
            ymin = -0.05
            ymax = 0.1
        interp.csvcompexp(sys.argv[1:], 'Flow around a conical island, case ' + scen[1], 'csvcomp', gauges, \
            'gauges_' + scen[1] + '.csv', 1.0, -0.32, -20.0, 6.0, 20.0, ymin, ymax)
