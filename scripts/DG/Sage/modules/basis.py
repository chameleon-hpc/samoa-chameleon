import numpy as np
from sage.functions.other import imag
from sage.calculus.var import var
from sage.rings.real_mpfr import RR
from sage.calculus.functional import simplify,derivative
from sage.matrix.constructor import matrix,vector,Matrix
from sage.functions.orthogonal_polys import gen_legendre_P

#from: http://mathworld.wolfram.com/LobattoQuadrature.html
def getGLLNodes(order):
    t = var("t",domain=RR)
    nodes = []
    if order > 1:
        print(derivative(gen_legendre_P(order,0,t),t))
        roots = np.roots(derivative(gen_legendre_P(order,0,t),t).coefficients(sparse=False)[::-1])
        #roots = simplify(derivative(gen_legendre_P(order,0,t),t).roots())
        nodes = [ node for node in roots if imag(node) == 0]
    nodes.extend([-1,1])
    nodes = sorted([ node/2.0 + 0.5 for node in nodes])
    #nodes = [node.numerical_approx() for node in nodes]
    return nodes

def getAlphaNodes(order):
    with open("alpha_nodes/nodes_"+str(order)+".txt") as nodes_file:
        nodes = nodes_file.read().split("\n")
        nodes = [node.split(" ") for node in nodes]
        nodes_float = []
        for node in nodes:
            node_float = []
            append=True
            for coord in node:
                try:
                    node_float.append(float(coord))
                except ValueError:
                    append=False
                    print("Ignoring :"+str(node))
            if append:
                nodes_float.append(node_float)
    return nodes_float

def generateAlphaWeights(basis):
    weights = [0.0]*len(basis)
    for i in range(len(basis)):
        weights[i] = integrate(integrate(basis[i],y,a=0,b=1.0-x),x,a=0.0,b=1.0)
    return weights
    

def getLagrangeBasis(nodes):
    basisTime = [0] * len(nodes)
    R = RR["t"]
    for i in range(0,len(nodes)):
        vals = [0] * len(nodes)
        vals[i] = 1
        tuples = zip(nodes, vals)

        basisTime[i] = R.lagrange_polynomial(tuples)
    return basisTime

#from: http://mathworld.wolfram.com/LobattoQuadrature.html
def getGLLNodes(order):
    t = var("t",domain=RR)
    nodes = []
    if order > 1:
        print(derivative(gen_legendre_P(order,0,t),t))
        roots = np.roots(derivative(gen_legendre_P(order,0,t),t).coefficients(sparse=False)[::-1])
        #roots = simplify(derivative(gen_legendre_P(order,0,t),t).roots())
        nodes = [ node for node in roots if imag(node) == 0]
    nodes.extend([-1,1])
    nodes = sorted([ node/2.0 + 0.5 for node in nodes])
    #nodes = [node.numerical_approx() for node in nodes]
    return nodes
