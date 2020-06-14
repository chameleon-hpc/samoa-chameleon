# Samoa^2 Alpha optimized nodes helper

# This module reads the files from scripts/DG/Sage/alpha_nodes and provides them as structured objects

import os, math
import pandas as pd


class NodeInformation:
    """This class holds values for a single node"""

    def __init__(self, x, y, w):
        self.x = x
        self.y = y
        self.w = w

    def __str__(self):
        return "X: " + "{0:0.20f}".format(self.x) + ", Y: " + \
            "{0:0.20f}".format(self.y) + ", Weight: " + "{0:0.20f}".format(self.w)

class AlphaNodes:
    """This class holds the alpha optimized nodes location and weights"""

    def __init__(self):
        basepath = os.path.dirname(__file__) + "/../../DG/Sage/alpha_nodes/"
        self.configs = [None] * 10
        self.configs[0] = [ NodeInformation(0.0, 0.0, 0.5) ]
        for n in range(1, 10):
            df = pd.read_csv(basepath + 'nodes_' + str(n) + '.txt', sep=' ', header=None)
            df.columns = ['x', 'y']
            dfw = pd.read_csv(basepath + 'weights_' + str(n) + '.txt', sep=' ', header=None)
            dfw.columns = ['w']
            df['w'] = dfw['w']
            entries = []
            for _, row in df.iterrows(): entries.append(NodeInformation(row['x'], row['y'], row['w']))
            self.configs[n] = entries

    def __str__(self):
        retstr = "Alpha optimized nodes information:\n"
        for i, config in enumerate(self.configs):
            retstr += " - Order: " + str(i) + ", nodes: " + str(len(config))
            if i < len(self.configs) - 1:
                retstr += "\n"
        return retstr

    def order_to_dof(self, order):
        return int(0.5 * (order + 1) * (order + 2))

    def dof_to_order(self, dof):
        return int(0.5 * (-3.0 + math.sqrt((8.0 * dof) + 1.0)))
