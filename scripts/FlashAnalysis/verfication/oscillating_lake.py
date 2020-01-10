#!/usr/bin/python3

# Verification of the oscillating lake scenario in [Vater, 2019]

import sys, core
from signal import *
from numba import njit
from pysamoaxdmf.reader import Reader, Parameter
import matplotlib, re
matplotlib.use('svg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import itertools, glob
from os.path import basename, dirname
import interp, lerror, math
from time import time
from multiprocessing import Pool 

omega = math.sqrt(0.2 * 9.80616)
P = (2.0 * math.pi) / omega

# P =           4.485701465466374
# 2 * P =       8.971402930932747

ts = [2.0 * P]
zposl = [0]

xnew = np.linspace(-2.0, 2.0, 1000) # Interpolation points

levelsh = np.linspace(0.0, 0.1, 20) # Color levels for h
levelshu = np.linspace(-0.003, 0.0015, 20) # Color levels for hu
levelshv = np.linspace(0.0, 0.08, 20) # Color levels for hv


# Analytical ideal results
@njit
def analytic_b(x0, x1):
    return 0.1 * ((x0 ** 2) + (x1 ** 2))

@njit
def analytic_h(x0, x1, t):
    return max((0.1 * ( \
        (x0 * math.cos(omega * t)) + (x1 * math.sin(omega * t)) + 0.75 \
        )) - analytic_b(x0, x1), 0.0)

@njit
def analytic_hu(x0, x1, t):
    u = analytic_u(x0, x1, t)
    h = analytic_h(x0, x1, t)
    return (u[0] * h, u[1] * h)

@njit
def analytic_u(x0, x1, t):
    if analytic_h(x0, x1, t) > 0:
        factor = omega / 2.0
        return (-math.sin(omega * t) * factor, \
                math.cos(omega * t) * factor)
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
    def huy(t): return list(map(lambda x: (analytic_hu(x, 0, t))[1], xnew))
    def uy(t): return list(map(lambda x: (analytic_u(x, 0, t))[1], xnew))

    # Read data file and load steps
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    print(str(xdmf))
    stepidxs = list(map(xdmf.step_near_time, ts))
    for sx in zip(ts, stepidxs):
        xdmf.steps[sx[1]].load()
        print(str(xdmf.steps[sx[1]]))
        print('   index: ' + str(sx[1]) + ', time deviation: ' + str(xdmf.steps[sx[1]].time - sx[0]) + 's')
    for i, sx in enumerate(stepidxs):
        ts[0] = xdmf.steps[sx].time

    # Extract data arrays from reader and reshape them
    pts, bs, hs, hbs, huxs, uxs, huys, uys = interp.read_reshape(xdmf, stepidxs)

    # Interpolate across line over domain for all variables
    print('Interpolating...')
    results = list(map(lambda zpos: interp.interp((hbs, huxs, huys, uxs, uys, bs), zpos, pts, xnew), zposl))
    results = results[0]

    # Plot results
    plt.rcParams['figure.figsize'] = [20, 8]
    plt.suptitle(name + "\n" + 'Comparison to exact solution at z=' + str(zposl).strip('[]') + ', t=2P')
    print("Crosssection at z=" + str(zposl[0]))
    t = xdmf.steps[stepidxs[0]].time
    for i, namep in enumerate([("h(x) + b(x)", h, 1), ("hu(x)", hux, 2), ("hv(x)", huy, 3), ("u(x)", ux, 5), ("v(x)", uy, 6)]):
        ax = plt.subplot(2, 3, namep[2])
        ax.set_title(namep[0] + ' at t=' + str(t) + ', z=' + str(zposl[0]))
        plt.axhline(0, color='grey', linestyle=':')
        plt.axvline(0, color='grey', linestyle=':')
        if i == 0: plt.plot(xnew, results[5][0], linestyle='--', color='g')
        plt.plot(xnew, namep[1](t), color='r')
        plt.plot(xnew, results[i][0], color='b')
        plt.xlim(-2, 2)
        if i == 0: plt.ylim(0, 0.4)
        elif i == 1: plt.ylim(-0.004, 0.0015)
        elif i == 2: plt.ylim(-0.005, 0.08)
        elif i == 3: plt.ylim(-0.05, 0.15)
        elif i == 4: plt.ylim(-0.05, 0.8)

    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(16, 8)
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

    t = ts[0] # Select t = 2P
    sx = xdmf.step_near_time(t)
    xdmf.steps[sx].load()
    print(str(xdmf.steps[sx]))
    print('   index: ' + str(sx) + ', time deviation: ' + str(xdmf.steps[sx].time - t) + 's')
    t = xdmf.steps[sx].time

    # Extract data arrays from reader and reshape them
    pts, bs, hs, hbs, huxs, _, huys, _ = interp.read_reshape(xdmf, [sx])

    # Filter dry/wet
    data = list(filter(lambda pack: pack[1] > dw, zip(pts[0], hs[0], huxs[0], huys[0])))
    ptsf, hsf, huxsf, huysf = zip(*(data))
    xs, zs = zip(*(ptsf))

    # Plot results
    plt.rcParams['figure.figsize'] = [17, 5]
    fig, (axh, axhu, axhv) = plt.subplots(1, 3)
    fig.suptitle(name + "\n" + 'Top-down contour plots at t=2P, dry/wet tolerance=' + str(float(scenario)))

    axh.set_title('h(x) at t=' + str(t))
    axh.set(xlim=(-2, 2), ylim=(-2, 2))
    axh.set_aspect('equal')
    axh.tricontour(xs, zs, hsf, levels=levelsh, linewidths=0.5, colors='k', extend='both')
    cntrh = axh.tricontourf(xs, zs, hsf, levels=levelsh, cmap='jet', extend='both')
    fig.colorbar(cntrh, ax=axh, shrink=0.5)

    axhu.set_title('hu(x) at t=' + str(t))
    axhu.set(xlim=(-2, 2), ylim=(-2, 2))
    axhu.set_aspect('equal')
    axhu.tricontour(xs, zs, huxsf, levels=levelshu, linewidths=0.5, colors='k', extend='both')
    cntrhu = axhu.tricontourf(xs, zs, huxsf, levels=levelshu, cmap='jet', extend='both')
    fig.colorbar(cntrhu, ax=axhu, shrink=0.5)

    axhv.set_title('hv(x) at t=' + str(t))
    axhv.set(xlim=(-2, 2), ylim=(-2, 2))
    axhv.set_aspect('equal')
    axhv.tricontour(xs, zs, huysf, levels=levelshv, linewidths=0.5, colors='k', extend='both')
    cntrhv = axhv.tricontourf(xs, zs, huysf, levels=levelshv, cmap='jet', extend='both')
    fig.colorbar(cntrhv, ax=axhv, shrink=0.5)

    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(16, 5)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')

def anacompseries(argv, name, limiter, filename):
    # Use getopt to parse argv
    inputdir, _, outputdir, _, _ = core.args(argv)
    print('Reading XMF files from dir: ' + inputdir)
    print('Writing output analytical comparison graphs to: ' + outputdir)
   
    seriesdirs = sorted(glob.glob(inputdir + '/oscillating_lake-*-light_' + limiter))
    p = re.compile('^.*oscillating_lake-([0-9]*)-light_BJ_.*$')
    seriesidxs = list(map(lambda x: int(p.findall(x)[0]), seriesdirs))
    deltax = list(map(lambda x: 4.0 / (2.0 ** (x / 2.0)), seriesidxs))

    # Aquire data
    l2_1d = []
    lsup_1d = []
    l2_2d = []
    lsup_2d = []
    dw = 0.00000001
    for sdir in seriesdirs:
        inputfile = glob.glob(sdir + '/*.xmf')[0]
        print('Reading XMF file: ' + inputfile)
        xdmf = Reader(inputfile, Parameter(2, 3, 3))
        print(str(xdmf))
        time = xdmf.steps[0].time
        #time = 0.0
        print('Loading step at t='+str(time))
        xdmf.steps[0].load()
        print('Calculating errors')
        l2_1d.append(math.sqrt(lerror.l2_1(xdmf.steps[0].cells_x, xdmf.steps[0].cells_h, time, analytic_h, dw)))
        lsup_1d.append(lerror.lsup_1(xdmf.steps[0].cells_x, xdmf.steps[0].cells_h, time, analytic_h, dw))
        l2_2d.append(math.sqrt(lerror.l2_2(xdmf.steps[0].cells_x, xdmf.steps[0].cells_h, xdmf.steps[0].cells_hu, time, analytic_hu, dw)))
        lsup_2d.append(lerror.lsup_2(xdmf.steps[0].cells_x, xdmf.steps[0].cells_h, xdmf.steps[0].cells_hu, time, analytic_hu, dw))
        xdmf.steps[0].unload()

    # Plot results
    plt.rcParams['figure.figsize'] = [14, 4.5]
    plt.suptitle(name + "\n" + r'$L^\infty$ and $L^2$ errors for $h$ and $hu$')
    plt.subplot(1, 2, 1)
    plt.plot(deltax, lsup_1d, marker='s')
    plt.plot(deltax, l2_1d, marker='o')
    plt.xlabel(r'$\Delta$x')
    plt.ylabel(r'$L^\infty$ (square) and $L^2$ (circle) errors $h$')
    plt.xscale('log')
    plt.yscale('log')
    plt.subplot(1, 2, 2)
    plt.plot(deltax, lsup_2d, marker='s')
    plt.plot(deltax, l2_2d, marker='o')
    plt.xlabel(r'$\Delta$x')
    plt.ylabel(r'$L^\infty$ (square) and $L^2$ (circle) errors error $hu$')
    plt.xscale('log')
    plt.yscale('log')

    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(14, 4.5)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')

def massenergy(argv, name, mode, filename):
    # Use getopt to parse argv
    inputfile, _, outputdir, _, workers = core.args(argv)
    print('Reading XMF file: ' + inputfile)
    print('Writing output ' + mode + ' graphs to: ' + outputdir)
    print('Using ' + str(workers) + ' worker(s)')

    # Read data file and load steps
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    pool = Pool(workers, init_proc_workers, [xdmf])
   
    # Main computation loop
    num_steps = len(xdmf.steps)
    t = [0.0] * num_steps
    data = [0.0] * num_steps
    starttime = time()
    proc = compute_proc_mass
    if(mode == "energy"): proc= compute_proc_energy
    for i, pdata in enumerate(pool.imap(proc, range(len(xdmf.steps))), 0):
        t[i], data[i] = pdata
        sys.stdout.write("\rComputing " + mode + "... n=" + str(i) + ", " + "{0:.2f}".format((float(i) / float(len(xdmf.steps))) * 100.0) + '%')
    endtime = time()
    print("\rComputed " + mode + " in " + "{0:.2f}".format(endtime - starttime) + 's            ')
    datanp = np.array(data)
    m0 = datanp[0]
    datanp = (datanp - m0) / m0

    # Plot results
    plt.rcParams['figure.figsize'] = [20, 8]
    plt.suptitle(name + "\n" + 'Total ' + mode + ' and relative ' + mode + ' change')
    plt.subplot(1, 2, 1)
    if(mode == "mass"): plt.ylabel("total mass [kg]")
    else: plt.ylabel("total energy [J]")
    plt.xlabel("t [s]")
    plt.plot(t, data)
    plt.subplot(1, 2, 2)
    if(mode == "mass"): plt.ylabel("Relative mass change [kg] (M - M0)/M0")
    else: plt.ylabel("relative energy change [J] (E - E0)/E0")
    plt.xlabel("t [s]")
    plt.plot(t, datanp)
    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(20, 8)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')

xdmf_p = None
def init_proc_workers(xdmf):
    global xdmf_p
    # Initialize worker process env
    xdmf_p = xdmf

def compute_proc_mass(i):
    global xdmf_p
    # Load a single step and compute its mass
    xdmf_p.steps[i].load()
    result = (xdmf_p.steps[i].time, lerror.mass(xdmf_p.steps[i].cells_x, xdmf_p.steps[i].cells_h, 0.001))
    xdmf_p.steps[i].unload()
    return result

def compute_proc_energy(i):
    global xdmf_p
    # Load a single step and compute its energy
    xdmf_p.steps[i].load()
    result = (xdmf_p.steps[i].time, lerror.energy(xdmf_p.steps[i].cells_x, xdmf_p.steps[i].cells_h, \
        xdmf_p.steps[i].cells_hu, xdmf_p.steps[i].cells_b, 0.001))
    xdmf_p.steps[i].unload()
    return result


if __name__ == '__main__':
    signal(SIGINT, core.siginth)
    _, _, _, scenario, _ = core.args(sys.argv[1:])
    if scenario == "series_BJ_vertex":
        anacompseries(sys.argv[1:], 'Oscillatory flow: errors in fluid depth and momentum', 'BJ_vertex', 'series')
    elif scenario == "series_BJ_edge":
        anacompseries(sys.argv[1:], 'Oscillatory flow: errors in fluid depth and momentum', 'BJ_edge', 'series')
    elif scenario == "mass":
        massenergy(sys.argv[1:], 'Oscillatory flow: Time series of mass', 'mass', 'mass')
    elif scenario == "energy":
        massenergy(sys.argv[1:], 'Oscillatory flow: Time series of energy', 'energy', 'energy')
    elif scenario == "":
        anacompzero(sys.argv[1:], 'Oscillatory flow in a parabolic bowl, P=' + str(P) + 's', 'anacomp')
    else:
        contour(sys.argv[1:], 'Oscillatory flow in a parabolic bowl, P=' + str(P) + 's', 'contour', float(scenario))