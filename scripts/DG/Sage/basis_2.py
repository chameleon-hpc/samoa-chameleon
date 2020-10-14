
def basis_0(x,y):
  return 2.0*(x - 1.0)**2 - 2.0*(-1.0*x + 1.0)*y - (-1.9999999999999996*x + 0.3999999999999999)*y + 2.0*y**2 + 1.0*x - 0.6000000000000001*y - 0.9999999999999998

def basis_1(x,y):
  return -4.0*(x - 1.0)**2 + (-4.0*x + 0.8)*y - 4.0*x - 0.8*y + 4.0

def basis_2(x,y):
  return 2.0*(x - 1.0)**2 + 3.0*x - 2.0

def basis_3(x,y):
  return (4.44089209850063e-16)*(x - 1.0)**2 + 4.0*(-1.0*x + 1.0)*y - 4.0*y**2 + (3.33066907387547e-16)*x - 3.33066907387547e-16

def basis_4(x,y):
  return (2.22044604925031e-16)*(x - 1.0)**2 + (4.0*x - 0.8)*y + (2.22044604925031e-16)*x + 0.8*y - 2.77555756156289e-16

def basis_5(x,y):
  return (1.11022302462516e-16)*(x - 1.0)**2 - 2.0*(-1.0*x + 1.0)*y + (-2.0*x + 0.4)*y + 2.0*y**2 + (1.11022302462516e-16)*x + 0.6000000000000001*y - 1.11022302462516e-16

basis=[basis_0,basis_1,basis_2,basis_3,basis_4,basis_5] 

def interpolate(x,y,Q):
    assert(len(Q) == len(basis))
    result = 0
    for i in range(len(Q)):
        result = result + basis[i](x,y) * Q[i]
    return result