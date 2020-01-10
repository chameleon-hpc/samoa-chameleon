#!/usr/bin/python3

# Verification of the longwave basin scenario in [Vater, 2019]

import sys, core
from signal import *
from numba import njit
from pysamoaxdmf.reader import Reader, Parameter
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy import interpolate
import itertools
from os.path import basename, dirname
import interp, lerror, math

A = ((2500.0 ** 4) - (2000.0 ** 4)) / \
    ((2500.0 ** 4) + (2000.0 ** 4))
omega = math.sqrt(8.0 * 9.80616) / 2500.0
P = (2.0 * math.pi) / omega

# P =           1773.129191803609
# 1.5 * P =     2659.6937877054135
# 1.75 * P =    3102.97608565631575
# 2 * P =       3546.258383607218
# 0.25 * P =    443.28229795090225

ts = [1.5 * P, 1.75 * P, 2.0 * P]
zposl = [0]

xnew = np.linspace(-4000, 4000, 1000) # Interpolation points
hulimits = [(-0.025, 0.025), (-0.8, 0.8), (-0.06, 0.06)] # Chart y limits hu
ulimits = [(-0.25, 0.25), (-2, 2), (-1.1, 1.1)] # Chart y limits u

levelsh = np.linspace(0, 1.5, 16) # Color levels for h
levelsu = np.linspace(-3, 3, 15) # Color levels for u


# Analytical ideal results
@njit
def analytic_b(x0, x1):
    return ((x0 ** 2) + (x1 ** 2)) / (2500 ** 2)

@njit
def analytic_h(x0, x1, t):
    return max(( \
            (math.sqrt(1.0 - (A ** 2))) / \
            (1.0           - (A * math.cos(omega * t))) \
        ) - ( \
            (((x0 ** 2) + (x1 ** 2)) * (1.0 - (A ** 2))) / \
            ((2500 ** 2)             * ((1.0 - (A * math.cos(omega * t))) ** 2)) \
        ), 0.0)

@njit
def analytic_hu(x0, x1, t):
    u = analytic_u(x0, x1, t)
    h = analytic_h(x0, x1, t)
    return (u[0] * h, u[1] * h)

@njit
def analytic_u(x0, x1, t):
    if analytic_h(x0, x1, t) > 0:
        factor = (omega * A * math.sin(omega * t)) / \
                 (2.0 * (1.0 - (A * math.cos(omega * t))))
        return (x0 * factor, x1 * factor)
    else:
        return (0.0, 0.0)

def anacompzero(argv, name, filename):
    # Use getopt to parse argv
    inputfile, _, outputdir, _, _ = core.args(argv)
    print('Reading XMF file: ' + inputfile)
    print('Writing output analytical comparison graphs to: ' + outputdir)
   
    # Compute analytical reference
    def h(t): return list(map(lambda x: analytic_h(x, 0, t) + analytic_b(x, 0), xnew))
    def hux(t): return list(map(lambda x: (analytic_hu(x, 0, t))[0], xnew))
    def ux(t): return list(map(lambda x: (analytic_u(x, 0, t))[0], xnew))

    # Read data file and load steps
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    print(str(xdmf))
    stepidxs = list(map(xdmf.step_near_time, ts))
    for sx in zip(ts, stepidxs):
        xdmf.steps[sx[1]].load()
        print(str(xdmf.steps[sx[1]]))
        print('   index: ' + str(sx[1]) + ', time deviation: ' + str(xdmf.steps[sx[1]].time - sx[0]) + 's')
    for i, sx in enumerate(stepidxs):
        ts[i] = xdmf.steps[sx].time

    # Extract data arrays from reader and reshape them
    pts, bs, hs, hbs, huxs, uxs, _, _ = interp.read_reshape(xdmf, stepidxs)

    # Interpolate across line over domain for all variables
    print('Interpolating...')
    results = list(map(lambda zpos: interp.interp((hbs, huxs, uxs, bs), zpos, pts, xnew), zposl))

    # Plot results
    plt.rcParams['figure.figsize'] = [20, 4]
    plt.suptitle(name + "\n" + 'Comparison to exact solution at z=' + str(zposl).strip('[]') + ', t=1.5P, 1.75P, 2P')
    numdia = len(zposl) * len(ts)
    diaidx = 0
    for q, result in enumerate(zip(results, zposl)):
        print("Crosssection at z=" + str(result[1]))
        for k, data in enumerate(zip(result[0][0], result[0][1], result[0][2], result[0][3], hulimits, ulimits, stepidxs)):
            t = xdmf.steps[data[6]].time
            for i, namep in enumerate([("h(x) + b(x)", h), ("hu(x)", hux), ("u(x)", ux)]):
                ax = plt.subplot(numdia, 3, diaidx + 1)
                ax.set_title(namep[0] + ' at t=' + str(t) + ', z=' + str(result[1]))
                plt.axhline(0, color='grey', linestyle=':')
                plt.axvline(0, color='grey', linestyle=':')
                if i == 0: plt.plot(xnew, data[3], linestyle='--', color='g')
                plt.plot(xnew, namep[1](t), color='r')
                plt.plot(xnew, data[i], color='b')
                plt.xlim(-4000, 4000)
                if i == 0: plt.ylim(0, 3)
                elif i == 1: plt.ylim(data[4][0], data[4][1])
                elif i == 2: plt.ylim(data[5][0], data[5][1])
                diaidx += 1
    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(16, 4 * numdia)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')

def contour(argv, name, filename, dw):
    # Use getopt to parse argv
    inputfile, _, outputdir, scenario, _ = core.args(argv)
    print('Reading XMF file: ' + inputfile)
    print('Writing output contour graphs to: ' + outputdir)
   
    # Read data file and load steps
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    print(str(xdmf) + "\n")

    t = ts[2] # Select t = 2P
    sx = xdmf.step_near_time(t)
    xdmf.steps[sx].load()
    print(str(xdmf.steps[sx]))
    print('   index: ' + str(sx) + ', time deviation: ' + str(xdmf.steps[sx].time - t) + 's')
    t = xdmf.steps[sx].time

    # Extract data arrays from reader and reshape them
    pts, bs, hs, hbs, huxs, uxs, _, _ = interp.read_reshape(xdmf, [sx])

    # Filter dry/wet
    xs, zs = zip(*(pts[0]))
    triag = tri.Triangulation(xs, zs)
    triag.set_mask([np.any(hs[0][tri] < dw) for tri in triag.triangles])

    # Plot results
    plt.rcParams['figure.figsize'] = [16, 5]
    fig, (axh, axhu, axu) = plt.subplots(1, 3)
    fig.suptitle(name + "\n" + 'Top-down contour plots at t=2P, dry/wet tolerance=' + str(float(scenario)))
    bathch = plt.Circle((0, 0), 2000, color='fuchsia', fill=False, lw=2, ls='--')
    bathchu = plt.Circle((0, 0), 2000, color='fuchsia', fill=False, lw=2, ls='--')
    bathcu = plt.Circle((0, 0), 2000, color='fuchsia', fill=False, lw=2, ls='--')

    axh.set_title('h(x) at t=' + str(t))
    axh.set(xlim=(-4000, 4000), ylim=(-4000, 4000))
    axh.set_aspect('equal')
    axh.tricontour(triag, hs[0], levels=levelsh, linewidths=0.5, colors='k', extend='max')
    cntrh = axh.tricontourf(triag, hs[0], levels=levelsh, cmap='jet', extend='max')
    axh.add_artist(bathch)
    fig.colorbar(cntrh, ax=axh, shrink=0.5)

    axhu.set_title('hu(x) at t=' + str(t))
    axhu.set(xlim=(-4000, 4000), ylim=(-4000, 4000))
    axhu.set_aspect('equal')
    axhu.tricontour(triag, huxs[0], levels=10, linewidths=0.5, colors='k', extend='both')
    cntrhu = axhu.tricontourf(triag, huxs[0], levels=10, cmap='jet', extend='both')
    axhu.add_artist(bathchu)
    fig.colorbar(cntrhu, ax=axhu, shrink=0.5)

    axu.set_title('u(x) at t=' + str(t))
    axu.set(xlim=(-4000, 4000), ylim=(-4000, 4000))
    axu.set_aspect('equal')
    axu.tricontour(triag, uxs[0], levels=levelsu, linewidths=0.5, colors='k', extend='both')
    cntru = axu.tricontourf(triag, uxs[0], levels=levelsu, cmap='jet', extend='both')
    axu.add_artist(bathcu)
    fig.colorbar(cntru, ax=axu, shrink=0.5)

    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(16, 5)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')

if __name__ == '__main__':
    signal(SIGINT, core.siginth)
    _, _, _, scenario, _ = core.args(sys.argv[1:])
    if scenario == "":
        anacompzero(sys.argv[1:], 'Long wave resonance in a paraboloid basin, P=' + str(P) + 's', 'anacomp')
    else:
        contour(sys.argv[1:], 'Long wave resonance in a paraboloid basin, P=' + str(P) + 's', 'contour', float(scenario))