from numba import njit

@njit
def basis_0(x,y):
  return -1.0*x - 1.0*y + 1.0

@njit
def basis_1(x,y):
  return 1.0*x - 5.55111512312578e-17

@njit
def basis_2(x,y):
  return 1.0*y - 5.55111512312578e-17

@njit
def interpolate(x, Q):
	assert(len(Q) == 3)
	return \
		(basis_0(x[0], x[1]) * Q[0]) + \
		(basis_1(x[0], x[1]) * Q[1]) + \
		(basis_2(x[0], x[1]) * Q[2])
