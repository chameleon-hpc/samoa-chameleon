import os
from SCons.Variables import *

def generate (env, **kw):

    SWE_dir="src/SWE/"
    
    matrices_template=SWE_dir+"SWE_dg_matrices.template"
    matrices_out=SWE_dir+"SWE_dg_matrices.f90"

    with open(matrices_template, "rt") as fin:
        with open(matrices_out, "wt") as fout:
            for line in fin:
                fout.write(line.replace('SWE_DG_ORDER_TAG', str(env['swe_dg_order'])))

    
    print("generate Kernels Tool")
    return

def exists(env):
    return True
