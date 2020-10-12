import numpy as np
from numba import njit
import math
from .basis import basis_1, basis_2, basis_3, basis_4, basis_5, basis_6, basis_7, basis_8


@njit
def _order_to_dof(order):
    return int(0.5 * (order + 1) * (order + 2))

@njit
def _dof_to_order(dof):
    return int(0.5 * (-3.0 + math.sqrt((8.0 * dof) + 1.0)))

@njit
def tri_cart_to_bary(x, v):
    det = ((v[1][1] - v[2][1]) * (v[0][0] - v[2][0])) + \
            ((v[2][0] - v[1][0]) * (v[0][1] - v[2][1]))
    l1 = (((v[1][1] - v[2][1]) * (x[0] - v[2][0])) + \
            ((v[2][0] - v[1][0]) * (x[1] - v[2][1]))) / det
    l2 = (((v[2][1] - v[0][1]) * (x[0] - v[2][0])) + \
            ((v[0][0] - v[2][0]) * (x[1] - v[2][1]))) / det
    return (l1, l2, 1.0 - l1 - l2)


@njit
def _rotate(x, t):
    return (np.cos(t) * x[0] - np.sin(t) * x[1], 
        np.sin(t) * x[0] + np.cos(t) * x[1])

@njit
def tri_cart_to_uni(x, v, plotter):
    vt = v - v[0]
    angle = -(math.pi / 4.0 * abs(plotter))
    # print(180.0 * angle / math.pi)

    for i in range(0, 3): vt[i] = _rotate(vt[i], angle)
    dist = max(abs(vt[0][0] - vt[1][0]), abs(vt[0][1] - vt[1][1]))
    # for i in range(0, 3): vt[i] = np.round(vt[i] / dist, 8)

    # print("dist=" + str(dist))

    # print("after rotate:")
    # print(vt)
    xo = (x[0] - v[0][0], x[1] - v[0][1])
    xn = _rotate(xo, angle)

    if(plotter > 0): return (xn[1] / dist, xn[0] / dist)
    else: return (xn[0] / dist, xn[1] / dist)


@njit
def tri_contains(x, v):
    res = tri_cart_to_bary(x, v)
    return res[0] >= 0 and res[1] >= 0 and res[2] >= 0

@njit
def tri_interp_bary_linear(bx, q):
    return (bx[0] * q[0]) + (bx[1] * q[1]) + (bx[2] * q[2])

# @njit
def tri_interp_cartuni_lagrange(bx, q):
    order = _dof_to_order(len(q))
    if order == 1:
        return basis_1.interpolate(bx[0], bx[1], q)
    elif order == 2:
        return basis_2.interpolate(bx[0], bx[1], q)
    elif order == 3:
        return basis_3.interpolate(bx[0], bx[1], q)
    elif order == 4:
        return basis_4.interpolate(bx[0], bx[1], q)
    elif order == 5:
        return basis_5.interpolate(bx[0], bx[1], q)
    elif order == 6:
        return basis_6.interpolate(bx[0], bx[1], q)
    elif order == 7:
        return basis_7.interpolate(bx[0], bx[1], q)
    elif order == 8:
        return basis_8.interpolate(bx[0], bx[1], q)
    else:
        return None

