# Samoa^2 Alpha exact optimized nodes helper

# This module reads the files from ./exact and provides them as structured objects

import os, math
import pandas as pd
import numpy as np

class ExactQuadRule:
    """This class holds the exact alpha optimized nodes location and weights"""

    def __init__(self):
        basepath = os.path.dirname(__file__) + "/exact/"
        self.configs = [None] * 10
        self.transformations = [None] * 10
        self.configs[0] = [ [1.0 / 3.0, 1.0 / 3.0, 0.5] ]
        for n in range(2, 8):
            self.configs[n]        = np.load(basepath+"/quad_"+str(n)+".npy")
            self.configs[n][:,2]   = self.configs[n][:,2] * 0.5
            self.transformation[n] = np.load(basepath+"/transform_"+str(n)+".npy")

    def __str__(self):
        retstr = "Exact Quadrature Rule Information:\n"
        for i, config in enumerate(self.configs):
            retstr += " - Order: " + str(i) + ", nodes: " + str(len(config))
            if i < len(self.configs) - 1:
                retstr += "\n"
        return retstr

    def order_to_dof(self, order):
        return int(0.5 * (order + 1) * (order + 2))

    def dof_to_order(self, dof):
        return int(0.5 * (-3.0 + math.sqrt((8.0 * dof) + 1.0)))
<<<<<<< HEAD
=======

>>>>>>> a6fd7dcce4fc601164dad7e3eeaedf752a16d985
