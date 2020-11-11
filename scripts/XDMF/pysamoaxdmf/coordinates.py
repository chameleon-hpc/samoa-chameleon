#!/usr/bin/env python3

# Coordinate transformation functions for FLASH and DG cells

import numpy as np
from numba import njit

@njit
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
def _tri_cart_to_uni_t(x, v):
    """This function computes a transformation matrix to convert a cartesian point on a DG cell to a normalized cartesian one on the unit triangle"""
    # x' = T*x with T = I_trig*(V^-1)
    uni = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]])
    tri = np.array([[v[0][0], v[1][0], v[2][0]], [v[0][1], v[1][1], v[2][1]], [1.0, 1.0, 1.0]])
    return np.dot(uni, np.linalg.inv(tri))

@njit
def tri_cart_to_uni(x, v):
    """This function converts a cartesian point on a DG cell to a normalized cartesian one on the unit triangle"""
    return np.dot(_tri_cart_to_uni_t(x, v), np.array([x[0], x[1], 1]))[0:2]

@njit
def tri_uni_to_cart(x, v):
    """This function converts a normalized cartesian point on the unit triangle to a cartesion one on a DG cell"""
    return np.dot(np.linalg.inv(_tri_cart_to_uni_t(x, v)), np.array([x[0], x[1], 1]))[0:2]

@njit
def tri_contains(x, v):
    """This function checks if a triangle contains a point"""
    res = tri_cart_to_bary(x, v)
    return res[0] >= 0 and res[1] >= 0 and res[2] >= 0