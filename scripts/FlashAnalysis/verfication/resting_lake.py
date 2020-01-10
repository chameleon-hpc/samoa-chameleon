#!/usr/bin/python3

# Verification of the resting lake scenario in [Vater, 2019]

import sys, core
from signal import *
import lerror
from numba import njit

# Analytical ideal results
@njit
def analytic_b_island(x0, x1):
    return max(0.0, 0.25 - (5.0 * \
        ( ((x0 - 0.5) ** 2) + ((x1 - 0.5) ** 2) )))

@njit
def analytic_h_island(x0, x1, t):
    return max(0.0, 0.1 - analytic_b_island(x0, x1))

@njit
def analytic_b_overlapping(x0, x1):
    om_1 = ((x0 - 0.35) ** 2) + ((x1 - 0.65) ** 2) < 0.01
    om_2 = ((x0 - 0.55) ** 2) + ((x1 - 0.45) ** 2) < 0.01
    om_3 = (abs(x0 - 0.47) < 0.25) and (abs(x1 - 0.55) < 0.25)
    om_4 = ((x0 - 0.5) ** 2) + ((x1 - 0.5) ** 2) < (0.45 ** 2)
    if om_1:
        return 0.15
    elif om_2:
        return 0.05
    elif om_3 and not(om_1 or om_2):
        return 0.07
    elif om_4 and not(om_3):
        return 0.03
    return 0.0

@njit
def analytic_h_overlapping(x0, x1, t):
    return max(0.0, 0.1 - analytic_b_overlapping(x0, x1))

@njit
def analytic_hu(x0, x1, t):
    return (0.0, 0.0)

if __name__ == '__main__':
    signal(SIGINT, core.siginth)
    _, _, _, scenario, _ = core.args(sys.argv[1:])
    if scenario == "island":
        lerror.anaerrors(sys.argv[1:], analytic_h_island, analytic_hu, 0.000001, 'Lake at Rest: First bathymetry (center island)', scenario)
    elif scenario == "overlapping":
        lerror.anaerrors(sys.argv[1:], analytic_h_overlapping, analytic_hu, 0.000001, 'Lake at Rest: Second bathymetry (overlapping)', scenario)
    else:
        print("Unknown scenario " + scenario)