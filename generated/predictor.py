#!/usr/bin/env python3
from yateto.gemm_configuration import *
from yateto.input import parseJSONMatrixFile
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
        DTDX = 1
        D_INDEX=2
        B_INDEX=2
        J_S_INDEX=SWE_DG_DOFS
        L_S_INDEX=SWE_DG_DOFS
        J_T_INDEX=SWE_DG_ORDER+1
        L_T_INDEX=SWE_DG_ORDER
        M_INDEX=3
        Q_INDEX=3
        l_s="l"
        j_s="j"
        l_t="L"
        j_t="J"
        d="d"
        b="b"
        m="m"
        q="q"
        
        max_depth=8
        # results
        P     = Tensor('P', (L_S_INDEX, L_T_INDEX, Q_INDEX))
        
        # input
        F     = Tensor('F',  ( J_S_INDEX, J_T_INDEX, B_INDEX, Q_INDEX ))
        S1    = Tensor('S1', ( J_S_INDEX, J_T_INDEX, M_INDEX ))
        S_res = Tensor('Sres', ( J_S_INDEX, J_T_INDEX, M_INDEX ))
        S2    = Tensor('S2', ( J_S_INDEX, J_T_INDEX, M_INDEX, D_INDEX))
        #J     = [Tensor('J_{}'.format(depth)    , ( Q_INDEX, M_INDEX)) for depth in range(1,max_depth+1)]
        clones={}
        db = parseJSONMatrixFile('{}/J_.json'.format(self.matrixDir),
                                 clones, alignStride=(lambda name: True))
        J_inv = [Tensor('J_inv_{}'.format(depth), (D_INDEX, D_INDEX)) for depth in range(1,max_depth+1)]
        Q0    = Tensor('Q0' ,( L_S_INDEX, Q_INDEX ))
        
        # precalculated matrices
        Kt11InvKt0               = Tensor('Kt11InvKt0'         , (L_T_INDEX,))
        Ktinv11_Mt               = Tensor('Ktinv11Mt'          , (J_T_INDEX,L_T_INDEX))
        Msinv_Ks                 = Tensor('MsinvKs'            , (J_S_INDEX,D_INDEX,L_S_INDEX))
        Msinv_Ks_Jinv            = Tensor('MsinvKs'            , (J_S_INDEX,B_INDEX,L_S_INDEX))    
        Ktinv11_Mt_J             = Tensor('Ktinv11MtJ'         , (J_T_INDEX,M_INDEX, Q_INDEX,L_T_INDEX))
        
        kernel_generator = lambda i : S_res[l_s+j_t+q] <= db.J[i][m+q] * S1[l_s+j_t+m]
        
        #    kernel_generator =  lambda i : P[l_s+l_t+q] <= P[l_s+l_t+q] + DTDX * (\
            #    Ktinv11_Mt[j_t+l_t]                       * db.J[i][m+q]   * S1 [l_s+j_t+m  ]+\
            #    Ktinv11_Mt[j_t+l_t] * Msinv_Ks[j_s+d+l_s] * db.J[i][m+q]   * S2 [j_s+j_t+m+d]+\
            #    Ktinv11_Mt[j_t+l_t] * Msinv_Ks_Jinv[j_s+b+l_s]       * F  [j_s+j_t+b+q])+\
            #    Kt11InvKt0[l_t] * Q0[l_s+q]
        
        g.addFamily('predictor', simpleParameterSpace(8),kernel_generator)
        
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

