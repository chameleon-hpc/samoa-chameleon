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


def interp(lists, zpos, pts, xnew):
    # Linear interpolation over a triangle mesh
    return list( \
        map(lambda data: \
            list(map(lambda pack: \
                interpolate.LinearNDInterpolator(pack[0], pack[1], 0) \
                    (list(zip(xnew, itertools.repeat(zpos)))), \
                zip(pts, data))), \
            lists))

def read_reshape(xdmf, stepidxs):
    # Extract data arrays from reader and reshape them
    pts = list(map(lambda sx: xdmf.steps[sx].cells_x.reshape(-1, xdmf.steps[sx].cells_x.shape[-1]), stepidxs))
    bs = list(map(lambda sx: xdmf.steps[sx].cells_b.reshape(-1), stepidxs))
    hs = list(map(lambda sx: xdmf.steps[sx].cells_h.reshape(-1), stepidxs))
    hbs = [a + b for a, b in zip(hs, bs)]
    huxs = list(map(lambda sx: (list(zip(*(xdmf.steps[sx].cells_hu.reshape(-1, xdmf.steps[sx].cells_hu.shape[-1])))))[0], stepidxs))
    huys = list(map(lambda sx: (list(zip(*(xdmf.steps[sx].cells_hu.reshape(-1, xdmf.steps[sx].cells_hu.shape[-1])))))[1], stepidxs))
    with np.errstate(divide='ignore', invalid='ignore'):
        uxs = [np.nan_to_num(np.divide(a, b)) for a, b in zip(huxs, hs)]
        uys = [np.nan_to_num(np.divide(a, b)) for a, b in zip(huys, hs)]
    return (pts, bs, hs, hbs, huxs, uxs, huys, uys)


def csvcompexp(argv, name, filename, gauges, gauges_filename, gauges_factor, gauges_shift, timeshift, xmin, xmax, ymin, ymax):
    # Use getopt to parse argv
    inputfile, datadir, outputdir, _, workers = core.args(argv)
    print('Reading XMF file: ' + inputfile)
    print('Reading CSV data from: ' + datadir)
    print('Writing output experimental (CSV) comparison graphs to: ' + outputdir)
    print('Using ' + str(workers) + ' worker(s)')
   
    # Load reference from csv files
    gauges_df = pd.read_csv(datadir + '/' + gauges_filename)
    gauges_pos_df = pd.read_csv(datadir + '/gauges_pos.csv')
    gauges_pos = {}
    for i in range(0, len(gauges)):
        g = gauges_pos_df.loc[gauges_pos_df['ch'] == gauges[i]]
        gauges_pos[gauges[i]] = np.array([g['x'].values[0], g['y'].values[0]])

    # Read data file and load steps
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    print(str(xdmf))
    gauge_data = {}
    for gauge in gauges: gauge_data[gauge] = [0.0] * len(xdmf.steps)
    time_data = [0.0] * len(xdmf.steps)
    for i in range(0, len(xdmf.steps)): time_data[i] = xdmf.steps[i].time

    # Interpolate across all time steps at gauge positions
    pool = Pool(workers, init_proc_interp, (xdmf, gauges, gauges_pos, gauges_shift))
    starttime = time()
    for i, pdata in enumerate(pool.imap(compute_proc_interp, range(len(xdmf.steps))), 0):
        for gauge in gauges: gauge_data[gauge][i] = pdata[gauge]
        sys.stdout.write("\rSampling gauge data... n=" + str(i) + ", " + "{0:.2f}".format((float(i) / float(len(xdmf.steps))) * 100.0) + '%')
    endtime = time()
    print("\rSampled gauge data in " + "{0:.2f}".format(endtime - starttime) + 's            ')

    # Plot results
    plt.rcParams['figure.figsize'] = [len(gauges) * 7, 7]
    fig, subplots = plt.subplots(1, len(gauges))
    fig.suptitle(name + "\n" + 'Comparison to experimental solution at gauges ' + str(gauges).strip('[]'))
    for i, subplot in enumerate(subplots):
        gauge = gauges[i]
        print('Plot for gauge #' + str(gauge))
        subplot.set_xlabel('Time [s]')
        subplot.set_ylabel('Water level [m]')
        subplot.set_title('Gauge ' + str(gauge) + ' at ' + str(gauges_pos[gauge]))
        subplot.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        subplot.plot(gauges_df["time"] + timeshift, gauges_df["ch" + str(gauge)] * gauges_factor, color='r')
        subplot.plot(time_data, gauge_data[gauge], color='b')

    print('Plotting...')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(len(gauges) * 7, 7)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')


xdmf_p = None
gauges_p = None
gauges_pos_p = None
gauges_shift_p = None
def init_proc_interp(xdmf, gauges, gauges_pos, gauges_shift):
    global xdmf_p, gauges_p, gauges_pos_p, gauges_shift_p
    # Initialize worker process env
    xdmf_p = xdmf
    gauges_p = gauges
    gauges_pos_p = gauges_pos
    gauges_shift_p = gauges_shift

def compute_proc_interp(j):
    global xdmf_p, gauges_p, gauges_pos_p, gauges_shift_p
    gauge_data_p = {}
    # Load a single step and compute gauge interpolation
    xdmf_p.steps[j].load()
    pts = xdmf_p.steps[j].cells_x.reshape(-1, xdmf_p.steps[j].cells_x.shape[-1])
    hs = xdmf_p.steps[j].cells_h.reshape(-1)
    bs = xdmf_p.steps[j].cells_b.reshape(-1)
    lvs = hs + bs
    interp = interpolate.LinearNDInterpolator(pts, lvs, 0.0)
    for gauge in gauges_p: gauge_data_p[gauge] = interp(gauges_pos_p[gauge])[0] + gauges_shift_p
    xdmf_p.steps[j].unload()
    return gauge_data_p