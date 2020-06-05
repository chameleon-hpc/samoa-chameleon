#!/usr/bin/env python3
from yateto.gemm_configuration import *
from yateto.input import parseJSONMatrixFile, parseJSONTensorFile
from yateto import *
from math import sqrt

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
        S     = Tensor('S' , (SWE_DG_DOFS, SWE_DG_ORDER+1, 3) )
        Q0    = Tensor('Q0', (SWE_DG_DOFS, 3) )
        
        clones={}
        db_m = parseJSONMatrixFile('{}/matrices_{}.json'.format(self.matrixDir,self.order),
                                   clones, alignStride=(lambda name: True))
        db_t = parseJSONTensorFile('{}/tensor_{}.json'.format(self.matrixDir,self.order),
                                   clones, alignStride=(lambda name: True))

        # precalculated matrices
        Msinv_Ks_Jinv = Tensor('MsinvKs', (SWE_DG_DOFS,2,SWE_DG_DOFS))    
        
        def predictor_generator(i):
            return P["lLq"] <= dtdx[""] *\
                ( db_m.t_k_t_11_inv_t_m_1["LJ"] * db_m.J[i]["mq"]                  * S ["lJm"]   +\
                  db_m.t_k_t_11_inv_t_m_1["LJ"] * db_t.s_m_inv_s_k_J_inv[i]["jbl"] * F ["jJbq"]) +\
                  db_m.t_k_t_11_inv_x_t_k_t_10["L"] * Q0["lq"]
         
#def kernel_generator(i):
#return P["lLq"] <= \
#                db_m.t_k_t_11_inv_t_m_1["LJ"]  * db_m.J[i]["mq"] * S ["lJm"]   +\
#                db_m.t_k_t_11_inv_x_t_k_t_10["L"] * Q0["lq"]

        
        
#        kernel_generator =  lambda i : return predictor_kernel(i)
        
        g.addFamily('predictor', simpleParameterSpace(8),predictor_generator)
        
        # kernel =  P[l_s+l_t+q] <= P[l_s+l_t+q] + DTDX * (\
            #                                                  Ktinv11_Mt[j_t+l_t]                       * J[m+q]   * S1 [l_s+j_t+m  ]+\
            #                                                  Ktinv11_Mt[j_t+l_t] * Msinv_Ks[j_s+d+l_s] * J[m+q]   * S2 [j_s+j_t+m+d]+\
            #                                                  Ktinv11_Mt[j_t+l_t] * Msinv_Ks_Jinv[j_s+b+l_s]       * F  [j_s+j_t+b+q])+\
            #                                                  Kt11InvKt0[l_t] * Q0[l_s+q]
        
        
        
        # kernel =  P[l_s+l_t+q] <= P[l_s+l_t+q] + DTDX * (\
            #                                                  Ktinv11_Mt[j_t+l_t]                       * J[m+q]   * S1 [l_s+j_t+m  ]+\
            #                                                  Ktinv11_Mt[j_t+l_t] * Msinv_Ks[j_s+d+l_s] * J[m+q]   * S2 [j_s+j_t+m+d]+\
            #                                                  Ktinv11_Mt[j_t+l_t] * Msinv_Ks_Jinv[j_s+b+l_s]       * F  [j_s+j_t+b+q])+\
            #                                                  Kt11InvKt0[l_t] * Q0[l_s+q]
        #                                                     Ktinv11_Mt[j_t+l_t]                       * J[q+m]     * S1 [l_s+j_t+m  ]+\
            #                                                     Ktinv11_Mt[j_t+l_t] * Msinv_Ks[j_s+d+l_s] * J[m+q]     * S2 [j_s+j_t+d+m]+\
            #                                                     Ktinv11_Mt_Msinv_Ks_Jinv[j_s+j_t+b+l_s+l_t]         * F  [j_s+j_t+q+b]+\    
        #g.add('predictor', kernel)

