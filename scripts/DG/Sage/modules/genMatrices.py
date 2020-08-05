from sage.calculus.var import var
from sage.rings.real_mpfr import RR
from sage.misc.functional import integrate
from sage.calculus.functional import simplify,derivative
from sage.functions.orthogonal_polys import jacobi_P
from sage.matrix.constructor import matrix,vector,Matrix
from sage.symbolic.constants import pi
from sage.functions.trig import cos,sin
from modules.meta import *

from sage.parallel.multiprocessing_sage import parallel_iter

from joblib import Parallel, delayed
import multiprocessing

def vandermonde_2D(order,nodes):   
    xi = var("xi",domain=RR)
    eta = var("eta",domain=RR)
 
    r = var("r",domain=RR)
    p = var("p",domain=RR)
    s = var("s",domain=RR)
    
    x = var("x",domain=RR)
    y = var("y",domain=RR)

    #span jacobi basis
    basis_list = []

    for d in range(0,order+1):
        for j in range(0,d+1):
            i = d-j
            poly = simplify(jacobi_P(i,0,0,r) * jacobi_P(j,2*i+1,0,s) *((1.0-s)*(1.0/2.0))**i)
            poly = poly( r = 2.0 * y / p - 1.0 ).simplify_full()
            poly = poly( s = -2.0*p+1.0 ).simplify_full()
            poly = poly( p = (1.0 - x) )
            basis_list.append(poly)
            
    basis= Matrix(basis_list[::-1])
    vandermonde = Matrix( [[0.0] * len(nodes)] * len(nodes))
    for i in range(0,len(nodes)):
        for j in range(0,len(nodes)):
            eta_n = nodes[i][0]
            xi_n  = nodes[i][1]
            
            vandermonde[i,j] = simplify(basis[0][j](x=eta_n,y=xi_n))
    vandermonde_inv = vandermonde.inverse()
    vandermonde_inv = [ [ entry if abs(entry) > 1.0e-15 else 0.0 for entry in row] for row in vandermonde_inv ]
    vandermonde_inv = matrix(vandermonde_inv)
    orthogonal_basis = simplify(basis * vandermonde_inv)
    return orthogonal_basis[0]

def generate_dudx(basis,direction):
    N = len(basis)
    dudx = Matrix([0]*N)
    for i in range(N):
        for j in range(N):
            dudx[i] = derivative(basis[i],direction)
    return dudx

def generate_s_m(basis):
    x = var("x",domain=RR)
    y = var("y",domain=RR)

    N = len(basis)
    m = Matrix([[0.0]*N]*N)
    for i in range(N):
        for j in range(N):
            #m[i,j] = integrate(integrate(simplify(basis[i]*basis[j]),y,a=0,b=1.0-x),x,a=0.0,b=1.0)
            #print(i)
            #print(j)
            #print(basis[i])
            m[i,j] = numerical_integral_2d(simplify(basis[i]*basis[j]))
    return m


def generate_t_m(basis):
    t = var("t",domain=RR)
    N = len(basis)
    m = Matrix([[0.0]*N]*N)
    for i in range(N):
        for j in range(N):
            m[i,j] = integrate((basis[i]*basis[j])(t),t,0,1)
    return m

def generate_t_a(basis):
    t = var("t",domain=RR)
    N = len(basis)
    m = Matrix([[0.0]*N])
    for i in range(N):
            m[0,i] = integrate(basis[i](t),t,0,1)
    return m

def generate_t_k(basis):
    N = len(basis)
    m = Matrix(RR,[[0.0]*N]*N)
    for i in range(N):
        for j in range(N):
            m[i,j] = integrate(basis[i]*derivative(basis[j],t)(t),t,0,1)
    return m

def generate_s_b(basis,order,edge,sparse=True):
    N = len(basis)
    if sparse:
        cols = getDGIndeces(order,edge)
    else:
        cols = range(0,N)
        
    b = Matrix([[0.0]*len(cols)]*N)
    t = var("t",domain=RR)
    if edge == "l":
        path = [t,0]
    if edge == "r":
        path = [0,t]
    if edge == "m":
        path = [1-t,t]
    for i in range(N):
        for j in range(len(cols)):
            b[i,j] = numerical_integral(basis[i](x=path[0],y=path[1]) * basis[cols[j]](x=path[0],y=path[1]), 0 , 1)[0]           
    return b

def generate_s_k(basis,direction):
    N = len(basis)
    k = Matrix([[0.0]*N]*N)
    for i in range(N):
        for j in range(N):
            k[i,j] = numerical_integral_2d(basis[j] * derivative(basis[i],direction))
    return k

def generate_phi(basis,order,triangles):
    N   = len(basis)
    phi = Matrix([[0.0]*N]*len(triangles))
    num=2*order+1
    for n in range(0,N):

        ind = 0
        for i in range(0,num):
            for j in range(0,i+1):
                if not j == 0:
                    ind = ind + 2
                triangle = triangles[ind]
                phi[ind,n] = integral(integral(basis[n],y,a=triangle[0][1],b=triangle[2][1]-(x-triangle[0][0])),
                                                        x,a=triangle[0][0],b=triangle[1][0])
            ind = ind + 1    
            
        ind = 2
        for i in range(1,num):
            for j in range(1,i+1):
                if not j == 1:
                    ind = ind + 2
                triangle = triangles[ind]
                            #triangles[ind] = [nodes_top[j], nodes_top[j-1], nodes_bottom[j]]
                phi[ind,n] = integral(integral(basis[n],x,a=triangle[0][0]-(y-triangle[2][1]),b=triangle[2][0]),
                                                        y,a=triangle[1][1],b=triangle[2][1])
            ind = ind + 3     

    return phi

def generate_mue(basis,order,phi):
    N = (order+1)*(order+2)//2
    N_mue = N + 1
    mue = matrix([[0.0]*N_mue]* N_mue)
    #print(phi.transpose()* phi)
    mue[:N,:N] = 2.0 * phi.transpose()* phi            
    for i in range(0,N):
        mue[N_mue-1,i] = numerical_integral_2d(basis[i])
        mue[i,N_mue-1] = mue[N_mue-1,i]
    return mue.inverse()


def printMatrix(matrix,name,order):
    dimString = dimMapping(matrix.nrows())+","+dimMapping(matrix.ncols())
    matrix_name = name+"("+dimString+")"
    output_string="real(kind=GRID_SR),Parameter :: "+matrix_name+" = reshape((/ &\n"
    for j in range(matrix.ncols()):
        for i in range(matrix.nrows()):
            #if i != 0 and i % 10 == 0:
            #        output_string = output_string + "&\n"
            output_string = output_string + str(matrix[i,j]) + "_GRID_SR ,"
        output_string = output_string + "&\n"
    #remove last seperator
    output_string = output_string[:-3]
    output_string = output_string + "  /),(/"+dimString+"/))"
    with open(str(name)+"_"+str(order)+".incl","w") as output:
        output.write(output_string)


def generate_ref(basis,side):
    #s_m = generate_s_m(basis)
    if(side == "l"):   #ref1
        trafo = [- 0.5 * x - 0.5 * y + 0.5,
                   0.5 * x - 0.5 * y + 0.5]
    elif(side == "r"): #ref2
        trafo = [- 0.5 * x + 0.5 * y + 0.5,
                 - 0.5 * x - 0.5 * y + 0.5]
    basis_trafo = [b(x=trafo[0], y=trafo[1]) for b in basis]
    N = len(basis)

    s_m_trafo = Matrix([[0.0]*N]*N)

    #old    
    #for i in range(N):
    #for j in range(N):
    #s_m_trafo[i,j] = numerical_integral_2d(simplify(basis[i]*basis_trafo[j]))

    for i in range(N):
        bases_arr = list(zip(range(N),[ simplify(basis[i]*basis_trafo[l]) for l in range(N)],[ gl_rule(N) for k in range(N)]))
        s_m_trafo_lst = sorted(list(numerical_integral_2d_arr(bases_arr)))
        s_m_trafo_i = map(lambda x : x[1], s_m_trafo_lst)
        s_m_trafo[i,:] = vector(s_m_trafo_i)

    return s_m_trafo

def generate_coarse(basis,side):
    #s_m = generate_s_m(basis)
    if(side == "l"):   #L
        trafo_a = [-(y - 0.5) + (x - 0.5),
                   -(y - 0.5) - (x - 0.5)]
        trafo_b = [y,x]
    elif(side == "r"): #R
        trafo_a = [-(x - 0.5) - (y - 0.5),
                   +(x - 0.5) - (y - 0.5)]
        trafo_b = [x,y]

    basis_trafo_a = [b(x=trafo_a[0], y=trafo_a[1]) for b in basis]
    basis_trafo_b = [b(x=trafo_b[0], y=trafo_b[1]) for b in basis]

    N = len(basis)
    threads = 8

    s_m_trafo_l = Matrix([[0.0]*N]*N)
    s_m_trafo_r = Matrix([[0.0]*N]*N)

    for i in range(N):
        bases_arr_l = list(zip(range(N),
                               [ simplify(basis_trafo_b[i]*basis_trafo_a[l]) for l in range(N)],
                               [ gl_rule(N,0.0,0.5,0.0,0.5,1) for k in range(N)]))
        
        s_m_trafo_lst = sorted(list(numerical_integral_2d_arr(bases_arr_l)))
        s_m_trafo_i = map(lambda x : x[1], s_m_trafo_lst)
        s_m_trafo_l[i,:] = vector(RR,s_m_trafo_i)


        bases_arr_r = list(zip(range(N),
                               [ simplify(basis_trafo_b[i]*basis_trafo_a[l]) for l in range(N)],
                               [ gl_rule(N,0.5,1.0,0.0,0.5,-1) for k in range(N)]))
                           
        s_m_trafo_lst = sorted(list(numerical_integral_2d_arr(bases_arr_r)))
        s_m_trafo_i = list(map(lambda x : x[1], s_m_trafo_lst))
        s_m_trafo_r[i,:] = vector(RR,s_m_trafo_i)

        #old
        #s_m_trafo[i,j+N] = numerical_integral_2d(simplify(basis[i]*basis_trafo[j]),0,0.5,0,0,1) + 
        # numerical_integral_2d(simplify(basis[i]*basis_trafo[j]),0.5,1,0,1,-1)

    return s_m_trafo_l + s_m_trafo_r



def getJacobian(i):
    angle = pi/4.0 * (i+1)
    return Matrix([[cos(angle),sin(angle)],[-sin(angle),cos(angle)]]).numerical_approx(digits=20)

def getJacobianPad(i):
    jac = getJacobian(i)
    return Matrix([[0,0],[jac[0,0],jac[1,0]],[jac[0,1],jac[1,1]]]).numerical_approx(digits=20)
