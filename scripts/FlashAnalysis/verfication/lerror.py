# L2 and Lsup error computation

from numba import njit
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
from time import time
from multiprocessing import Pool 
from functools import partial
from pysamoaxdmf.reader import Reader, Parameter
import math, sys, core
import math

@njit
def _tri_area(a, b, c):
    # Calculates the area of a 2D triangle
    ab = b - a
    ac = c - a
    return 0.5 * abs((ab[0] * ac[1]) - (ab[1] * ac[0]))

@njit
def _inner_1(num, ana, weight):
    return ((num - ana) ** 2) * weight

@njit
def _inner_0(num, ana):
    return abs(num - ana)

@njit
def _inner_2(num0, num1, ana0, ana1, weight):
    return (((num0 - ana0) 
    ** 2) + ((num1 - ana1) ** 2)) * weight

@njit
def l2_1(x, h, t, fana, dw):
    # Compute the L2 error for h
    cell_sum = 0.0
    for i in range(len(x)):
        inner_sum = 0.0
        weight = (1.0 / 3.0) * _tri_area(x[i][0], x[i][1], x[i][2])
        for j in range(3):
            h0 = 0.0 if (h[i][j] <= dw) else h[i][j]
            inner_sum += _inner_1(h0, fana(x[i][j][0], x[i][j][1], t), weight)
        cell_sum += inner_sum
    return cell_sum

@njit
def l2_2(x, h, hu, t, fana, dw):
    # Compute the L2 error for hu
    cell_sum = 0.0
    for i in range(len(x)):
        inner_sum = 0.0
        weight = (1.0 / 3.0) * _tri_area(x[i][0], x[i][1], x[i][2])
        for j in range(3):
            hu0 = 0.0 if (h[i][j] <= dw) else hu[i][j][0]
            hu1 = 0.0 if (h[i][j] <= dw) else hu[i][j][1]
            ana = fana(x[i][j][0], x[i][j][1], t)
            inner_sum += _inner_2(hu0, hu1, ana[0], ana[1], weight)
        cell_sum += inner_sum
    return cell_sum

@njit
def lsup_1(x, h, t, fana, dw):
    # Compute the Lsup error for h
    max_err = 0.0
    for i in range(len(x)):
        inner_max_err = 0.0
        for j in range(3):
            h0 = 0.0 if (h[i][j] <= dw) else h[i][j]
            inner_max_err = max(_inner_0(h0, fana(x[i][j][0], x[i][j][1], t)), inner_max_err)
        max_err = max(inner_max_err, max_err)
    return max_err

@njit
def lsup_2(x, h, hu, t, fana, dw):
    # Compute the Lsup error for hu
    max_err = 0.0
    for i in range(len(x)):
        inner_max_err = 0.0
        for j in range(3):
            hu0 = 0.0 if (h[i][j] <= dw) else hu[i][j][0]
            hu1 = 0.0 if (h[i][j] <= dw) else hu[i][j][1]
            ana = fana(x[i][j][0], x[i][j][1], t)
            hu_err = _inner_0(hu0, ana[0])
            hv_err = _inner_0(hu1, ana[1])
            inner_max_err = max(hu_err, inner_max_err)
            inner_max_err = max(hv_err, inner_max_err)
        max_err = max(inner_max_err, max_err)
    return max_err

@njit
def mass(x, h, dw):
    # Compute the total mass
    total_mass = 0.0
    for i in range(len(x)):
        cell_mass = 0.0
        weight = (1.0 / 3.0) * _tri_area(x[i][0], x[i][1], x[i][2])
        for j in range(3):
            h0 = 0.0 if (h[i][j] <= dw) else h[i][j]
            cell_mass += h0 * weight
        total_mass += cell_mass
    return total_mass

@njit
def energy(x, h, hu, b, dw):
    # Compute the total energy
    total_energy = 0.0
    for i in range(len(x)):
        cell_energy = 0.0
        weight = (1.0 / 3.0) * _tri_area(x[i][0], x[i][1], x[i][2])
        for j in range(3):
            h0 = 0.0 if (h[i][j] <= dw) else h[i][j]
            hu0 = 0.0 if (h[i][j] <= dw) else hu[i][j][0]
            hu1 = 0.0 if (h[i][j] <= dw) else hu[i][j][1]
            point_mass = h0 * weight
            cell_energy += point_mass * 9.80616 * ((h0 / 2.0) + b[i][j]) # potential energy
            cell_energy += 0.5 * point_mass * ((hu0 ** 2) + (hu1 ** 2)) # kinetic energy
        total_energy += cell_energy
    return total_energy


def anaerrors(argv, ana_h, ana_hu, dw, name, filename):
    # Use getopt to parse argv
    inputfile, _, outputdir, scenario, workers = core.args(argv)
    
    print('Reading XMF file: ' + inputfile)
    print('Writing output limiter error graphs to: ' + outputdir + '/' + filename)
    print('Using ' + str(workers) + ' worker(s)')
    xdmf = Reader(inputfile, Parameter(2, 3, 3))
    pool = Pool(workers, init_proc_lerror, (xdmf, ana_h, ana_hu, dw))

    # Main computation loop
    num_steps = len(xdmf.steps)
    t = [0.0] * num_steps
    l2e_h = [0.0] * num_steps
    lsupe_h = [0.0] * num_steps
    l2e_hu = [0.0] * num_steps
    lsupe_hu = [0.0] * num_steps
    starttime = time()
    for i, pdata in enumerate(pool.imap(compute_proc_lerror, range(len(xdmf.steps))), 0):
        t[i], l2e_h[i], lsupe_h[i], l2e_hu[i], lsupe_hu[i] = pdata
        sys.stdout.write("\rComputing limiter errors... n=" + str(i) + ", " + "{0:.2f}".format((float(i) / float(len(xdmf.steps))) * 100.0) + '%')
    endtime = time()
    print("\rComputed limiter errors in " + "{0:.2f}".format(endtime - starttime) + 's            ')

    # Plot data
    print('Plotting...')
    plt.suptitle(name + "\n" + r'$L^\infty$ and $L^2$ errors for $h$ and $hu$')
    plt.subplot(2, 2, 1)
    plt.plot(t, lsupe_h)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$L^\infty$ error $h$')
    plt.yscale('log')
    plt.subplot(2, 2, 2)
    plt.plot(t, l2e_h)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$L^2$ error $h$')
    plt.yscale('log')
    plt.subplot(2, 2, 3)
    plt.plot(t, lsupe_hu)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$L^\infty$ error $hu$')
    plt.yscale('log')
    plt.subplot(2, 2, 4)
    plt.plot(t, l2e_hu)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$L^2$ error $hu$')
    plt.yscale('log')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(12, 7.5)
    plt.savefig(outputdir + '/' + filename + '.svg', dpi=100, frameon=False)
    print('Done.')


xdmf_p = None
ana_h_p = None
ana_hu_p = None
dw_p = 0.0
def init_proc_lerror(xdmf, ana_h, ana_hu, dw):
    global xdmf_p, ana_h_p, ana_hu_p, dw_p
    # Initialize worker process env
    xdmf_p = xdmf
    ana_h_p = ana_h
    ana_hu_p = ana_hu
    dw_p = dw

def compute_proc_lerror(i):
    global xdmf_p, ana_h_p, ana_hu_p, dw_p
    # Load a single step and compute its errors
    xdmf_p.steps[i].load()
    result = (xdmf_p.steps[i].time, \
        math.sqrt(l2_1(xdmf_p.steps[i].cells_x, xdmf_p.steps[i].cells_h, xdmf_p.steps[i].time, ana_h_p, dw_p)), \
        lsup_1(xdmf_p.steps[i].cells_x, xdmf_p.steps[i].cells_h, xdmf_p.steps[i].time, ana_h_p, dw_p), \
        math.sqrt(l2_2(xdmf_p.steps[i].cells_x, xdmf_p.steps[i].cells_h, xdmf_p.steps[i].cells_hu, xdmf_p.steps[i].time, ana_hu_p, dw_p)), \
        lsup_2(xdmf_p.steps[i].cells_x, xdmf_p.steps[i].cells_h, xdmf_p.steps[i].cells_hu, xdmf_p.steps[i].time, ana_hu_p, dw_p))
    xdmf_p.steps[i].unload()
    return result
