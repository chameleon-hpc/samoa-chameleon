#!/usr/bin/env python3
from yateto.gemm_configuration import *
from yateto.input import parseJSONMatrixFile, parseJSONTensorFile
from yateto import *
from math import sqrt

class SourceGenerator():
    def __init__(self,order,matrixDir):
        self.order = int(order)
        self.matrixDir = matrixDir
    
    def gemm_cfg(self,arch):
        return GeneratorCollection([LIBXSMM(arch), MKL(arch)])


    def add(self,g):
        SWE_DG_ORDER = self.order
        SWE_DG_DOFS = (SWE_DG_ORDER+1)*(SWE_DG_ORDER+2)//2
        
        H  = Tensor('H',  (SWE_DG_DOFS, SWE_DG_ORDER+1, 2) )
        W  = Tensor('W',  (SWE_DG_DOFS, SWE_DG_ORDER+1) )
        Dx = Tensor('Dx', (SWE_DG_DOFS, SWE_DG_DOFS, 2) )
        S2 = Tensor('S2', (SWE_DG_DOFS, SWE_DG_ORDER+1, 2) )

        db_t = parseJSONTensorFile('{}/tensor_{}.json'.format(self.matrixDir,self.order),
                                   {}, alignStride=(lambda name: True), transpose=(lambda x: False))


        source = S2["lJm"] <= H["lJm"] * db_t.Dx["ljm"] * W["jJ"]

        g.add("compute_source",  source)    

class PredictorGenerator():
    def __init__(self,order,matrixDir):
        self.order = int(order)
        self.matrixDir = matrixDir
    
    def gemm_cfg(self,arch):
        return GeneratorCollection([LIBXSMM(arch), MKL(arch)])


    def add(self,g):
        SWE_DG_ORDER = self.order
        SWE_DG_DOFS = (SWE_DG_ORDER+1)*(SWE_DG_ORDER+2)//2
        max_depth=8

        dtdx  = Tensor("dtdx",())
        
        # results
        P     = Tensor('P', (SWE_DG_DOFS, SWE_DG_ORDER, 3))
        
        # input
        F     = Tensor('F' , (SWE_DG_DOFS, SWE_DG_ORDER+1, 2, 3) )
        S     = Tensor('S' , (SWE_DG_DOFS, SWE_DG_ORDER+1, 2) )
        Q0    = Tensor('Q0', (SWE_DG_DOFS, 3) )
       
        clones={}
        def transpose(name):
            dict_t = { "t_k_t_11_inv_t_m_1" : True}
            
            for mat in ["J","s_m_inv_s_k_J_inv"]:
                for n in ["{}({})".format(mat,i) for i in range(0,8)]:
                    dict_t[n] = True
            
            return dict_t[name] if name in dict_t else False

#        transpose=lambda x: False
        
        db_m = parseJSONMatrixFile('{}/matrices_{}.json'.format(self.matrixDir,self.order),
                                   clones, alignStride=(lambda name: True), transpose=transpose)
        db_t = parseJSONTensorFile('{}/tensor_{}.json'.format(self.matrixDir,self.order),
                                   clones, alignStride=(lambda name: True), transpose=transpose)

        def predictor_generator(i):
            return P["lLq"] <= dtdx[""] *\
                ( db_m.t_k_t_11_inv_t_m_1["JL"] * db_m.J[i]["mq"]                  * S ["lJm"]   +\
                  db_m.t_k_t_11_inv_t_m_1["JL"] * db_t.s_m_inv_s_k_J_inv[i]["jbl"] * F ["jJbq"]) -\
                db_m.t_k_t_11_inv_x_t_k_t_10["L"] * Q0["lq"]    

        g.addFamily('predictor', simpleParameterSpace(8),predictor_generator)

