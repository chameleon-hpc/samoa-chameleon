#!/usr/bin/env python3

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