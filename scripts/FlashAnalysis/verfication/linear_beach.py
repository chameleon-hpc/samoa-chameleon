#!/usr/bin/python3

# Verification of the linear beach scenario in [Vater, 2019]

import sys, core
from signal import *
from pysamoaxdmf.reader import Reader, Parameter
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import itertools
from os.path import basename, dirname
import pandas as pd
import interp


ts = [160, 175, 220]
zposl = [25000] #[24825, 24870, 24900] #[10, 25, 45] 

xnew = np.linspace(-400, 50000, 1000) # Interpolation points
hulimits = [(-90, 80), (-160, 10), (-270, 10)] # Chart y limits hu
ulimits = [(-1, 7), (-15, 11), (-4, 1)] # Chart y limits u


def csvcompzero(argv, name, filename):
    # Use getopt to parse argv
    inputfile, datadir, outputdir, _, _ = core.args(argv)
    print('Reading XMF file: ' + inputfile)
    print('Reading CSV data from: ' + datadir)
    print('Writing output analytical (CSV) comparison graphs to: ' + outputdir)
   
    # Load reference from csv files
    refdfs = list(map(lambda t: pd.read_csv(datadir + "/t" + str(t) + ".csv"), ts))
    for df in refdfs: df["hu"] = (df["eta"] + (0.1 * df["x"])) * df["u"]

    # Read data file and load steps
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    print(str(xdmf))
    stepidxs = list(map(xdmf.step_near_time, ts))
    for sx in zip(ts, stepidxs):
        xdmf.steps[sx[1]].load()
        print(str(xdmf.steps[sx[1]]))
        print('   index: ' + str(sx[1]) + ', time deviation: ' + str(xdmf.steps[sx[1]].time - sx[0]) + 's')
    for i, sx in enumerate(stepidxs):
        ts[i] = round(xdmf.steps[sx].time, 4)

    # Extract data arrays from reader and reshape them
    pts, bs, hs, hbs, huxs, uxs, _, _ = interp.read_reshape(xdmf, stepidxs)

    # Interpolate across line over domain for all variables
    print('Interpolating...')
    results = list(map(lambda zpos: interp.interp((hbs, huxs, uxs, bs), zpos, pts, xnew), zposl))

    # Plot results
    plt.rcParams['figure.figsize'] = [20, 4]
    plt.suptitle(name + "\n" + 'Comparison to exact solution at z=' + str(zposl).strip('[]') + ', t=' + str(ts).strip('[]'))
    numdia = len(zposl) * len(ts)
    diaidx = 0
    for q, result in enumerate(zip(results, zposl)):
        print("Crosssection at z=" + str(result[1]))
        for k, data in enumerate(zip(result[0][0], result[0][1], result[0][2], result[0][3], hulimits, ulimits, ts, refdfs)):
            for i, sub in enumerate(["eta", "hu", "u"]):
                ax = plt.subplot(numdia, 3, diaidx + 1)
                name = sub
                if i == 0: name = 'h(x) + b(x)'
                if i == 1: name = 'hu(x)'
                if i == 2: name = 'u(x)'
                ax.set_title(name + ' at t=' + str(data[6]) + ', z=' + str(result[1]))
                plt.axhline(0, color='grey', linestyle=':')
                plt.axvline(0, color='grey', linestyle=':')
                plt.plot(data[7]["x"], data[7][sub], color='r')
                plt.plot(xnew, data[i], color='b')
                if i == 0: plt.plot(xnew, data[3], linestyle='--', color='g')
                plt.xlim(-400, 800)
                if i == 0: plt.ylim(-25, 25)
                elif i == 1: plt.ylim(data[4][0], data[4][1])
                elif i == 2: plt.ylim(data[5][0], data[5][1])
                diaidx += 1
    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(16, 4 * numdia)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')


if __name__ == '__main__':
    signal(SIGINT, core.siginth)
    csvcompzero(sys.argv[1:], 'Tsunami runup onto a linearly sloping beach', 'csvcomp')
