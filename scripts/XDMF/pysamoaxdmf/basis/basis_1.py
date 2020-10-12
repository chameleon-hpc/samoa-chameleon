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

basis=[basis_0,basis_1,basis_2] 

def interpolate(x,y,Q):
    assert(len(Q) == len(basis))
    result = 0
    for i in range(len(Q)):
        result = result + basis[i](x,y) * Q[i]
    return result
