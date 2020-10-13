#!/usr/bin/python3

# Sampling functions for FLASH and DG cells

# This module provides functions for sampling values from a FLASH or DG cell
# and also some helper functions.
# It can also be invoked, and will then compile itself to a native library.

import math

import numpy as np
from numba import njit
from numba.pycc import CC

from pysamoaxdmf.basis import basis_1, basis_2, basis_3, basis_4, basis_5, basis_6, basis_7, basis_8

cc = CC('sampler')

@njit
def _order_to_dof(order):
    """This function calculates the amount of DG DOFs from the order"""
    return int(0.5 * (order + 1) * (order + 2))

@njit
def _dof_to_order(dof):
    """This function calculates the DG order from the amount of DOFs"""
    return int(0.5 * (-3.0 + math.sqrt((8.0 * dof) + 1.0)))

@njit
@cc.export('tri_cart_to_bary', '(f8[:], f8[:,:])')
def tri_cart_to_bary(x, v):
    """This function converts a cartesian point on a triangle to a barycentric one"""
    det = ((v[1][1] - v[2][1]) * (v[0][0] - v[2][0])) + \
            ((v[2][0] - v[1][0]) * (v[0][1] - v[2][1]))
    l1 = (((v[1][1] - v[2][1]) * (x[0] - v[2][0])) + \
            ((v[2][0] - v[1][0]) * (x[1] - v[2][1]))) / det
    l2 = (((v[2][1] - v[0][1]) * (x[0] - v[2][0])) + \
            ((v[0][0] - v[2][0]) * (x[1] - v[2][1]))) / det
    return np.array([l1, l2, 1.0 - l1 - l2])


@njit
def _rotate(x, t):
    """This function rotates a point around the origin using a specified angle"""
    return (np.cos(t) * x[0] - np.sin(t) * x[1], 
        np.sin(t) * x[0] + np.cos(t) * x[1])

@njit
@cc.export('tri_cart_to_uni', '(f8[:], f8[:,:], i8)')
def tri_cart_to_uni(x, v, plotter):
    """This function converts a cartesian point on a DG cell to a normalized one on the unit triangle"""
    vt = v - v[0]
    angle = -(math.pi / 4.0 * abs(plotter))
    for i in range(0, 3): vt[i] = _rotate(vt[i], angle)
    dist = max(abs(vt[0][0] - vt[1][0]), abs(vt[0][1] - vt[1][1]))
    xo = (x[0] - v[0][0], x[1] - v[0][1])
    xn = _rotate(xo, angle)
    if(plotter > 0): return np.array([xn[1] / dist, xn[0] / dist])
    else: return np.array([xn[0] / dist, xn[1] / dist])

@njit
@cc.export('tri_contains', '(f8[:], f8[:,:])')
def tri_contains(x, v):
    """This function checks if a triangle contains a point"""
    res = tri_cart_to_bary(x, v)
    return res[0] >= 0 and res[1] >= 0 and res[2] >= 0

@njit
def _tri_interp_bary_linear(bx, q):
    """This function interpolates a value on a FLASH cell"""
    return (bx[0] * q[0]) + (bx[1] * q[1]) + (bx[2] * q[2])

@njit
@cc.export('tri_interp_bary_linear_1d', '(f8[:], f8[:])')
def tri_interp_bary_linear_1d(bx, q):
    return _tri_interp_bary_linear(bx, q)

@njit
@cc.export('tri_interp_bary_linear_2d', '(f8[:], f8[:,:])')
def tri_interp_bary_linear_2d(bx, q):
    return _tri_interp_bary_linear(bx, q)


@njit
def _tri_interp_cartuni_lagrange(bx, q):
    """This function interpolates a value on a DG cell"""
    order = _dof_to_order(len(q))
    if order == 1:
        return basis_1.interpolate(bx, q)
    elif order == 2:
        return basis_2.interpolate(bx, q)
    elif order == 3:
        return basis_3.interpolate(bx, q)
    elif order == 4:
        return basis_4.interpolate(bx, q)
    elif order == 5:
        return basis_5.interpolate(bx, q)
    elif order == 6:
        return basis_6.interpolate(bx, q)
    elif order == 7:
        return basis_7.interpolate(bx, q)
    elif order == 8:
        return basis_8.interpolate(bx, q)
    else:
        return None

@njit
@cc.export('tri_interp_cartuni_lagrange_1d', '(f8[:], f8[:])')
def tri_interp_cartuni_lagrange_1d(bx, q):
    return _tri_interp_cartuni_lagrange(bx, q)

@njit
@cc.export('tri_interp_cartuni_lagrange_2d', '(f8[:], f8[:,:])')
def tri_interp_cartuni_lagrange_2d(bx, q):
    return _tri_interp_cartuni_lagrange(bx, q)


if __name__ == "__main__":
    print("Precompiling to native library, this might take a while.")
    cc.verbose = True
    cc.target_cpu = 'host'
    cc.compile()