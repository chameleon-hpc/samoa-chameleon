# Samoa^2 XDMF Reader for Python 3

# This module provides a easy-to-use object-orientated interface
# for the XDMF (HDF5 database and XML-based XMF index) data generated
# by the XDMF output module in samoa^2

import numpy as np
from numpy.testing.utils import HAS_REFCOUNT
import h5py
import os
import pyqtree
import sys
from lxml import etree
from . import sampler
import math


class Parameter:
    """This class represents the XDMF data and file layout configuration"""

    def __init__(self, valsg_width, valst_width, attr_width, attr_width_xml):
        self.valsg_width = valsg_width
        self.valst_width = valst_width
        self.attr_width = attr_width
        self.attr_width_xml = attr_width_xml

    def __str__(self):
        return str(self.valsg_width) + " - " + str(self.valst_width) + " - " + \
            str(self.attr_width) + " (" + str(self.attr_width_xml) + ")"

class StepLayer:
    """This class represens a sub-layer of a step, a set of cells with the same type (data layout)"""

    def __init__(self, param, name, h5path, h5gname, h5hname, h5huname, h5bname, h5lname):
        self.param = param
        self.name = name
        self.h5path = h5path
        self.h5gname = h5gname
        self.h5hname = h5hname
        self.h5huname = h5huname
        self.h5bname = h5bname
        self.h5lname = h5lname
        self.cells_length = 0
        self.unload()
        self.layer_name_map = {
            "c": "FV",
            "p": "DG"
        }

    def load(self):
        h5 = h5py.File(self.h5path, 'r')
        h5dg = h5[self.h5gname]
        h5dh = h5[self.h5hname]
        h5dhu = h5[self.h5huname]
        h5db = h5[self.h5bname]
        h5dl = h5[self.h5lname]
        # Actually load the data to speed up lookups
        gdata = h5dg[...]
        self.cells_length = h5dh.shape[0]
        self.cells_h = h5dh[...]
        self.cells_hu = h5dhu[...]
        self.cells_b = h5db[...]
        self.cells_l = h5dl[...]
        # Iterate all cells and unpack location
        self.cells_x = np.empty([self.cells_length, self.param.valst_width, gdata.shape[1]]) 
        for i in range(0, self.cells_length):
            self.cells_x[i] = gdata[(i * self.param.valst_width):(i * self.param.valst_width) + self.param.valst_width]
        h5.close()
        # Build quadtree index
        if self.cells_length > 0:
            bmin = np.min(gdata, 0)
            bmax = np.max(gdata, 0)
            self.index_bounds = (bmin[0], bmin[1], bmax[0], bmax[1])
            self.index = pyqtree.Index(bbox=self.index_bounds)
            for i, cell in enumerate(self.cells_x):
                cmin = np.min(cell, 0)
                cmax = np.max(cell, 0)
                self.index.insert(i, (cmin[0], cmin[1], cmax[0], cmax[1]))

    def unload(self):
        self.cells_length = 0
        self.cells_x = np.array([])
        self.cells_h = np.array([])
        self.cells_hu = np.array([])
        self.cells_b = np.array([])
        self.cells_l = np.array([])
        self.index = None
        self.index_bounds = (0, 0, 0, 0)

    def sample(self, x):
        if(self.index is None): 
            return ()
        # Lookup position in index
        candidates = self.index.intersect((x[0], x[1], x[0] + sys.float_info.epsilon, x[1] + sys.float_info.epsilon))
        if len(candidates) == 0:
            return ()
        for can in candidates:
            # Test candidate
            if(sampler.tri_contains(np.array(x), self.cells_x[can])):
                if(self.name == 'p'):
                    cuni = sampler.tri_cart_to_uni(np.array(x), self.cells_x[can], np.asscalar(self.cells_l[can]))
                    # DG cell sampling
                    return (sampler.tri_interp_cartuni_lagrange_1d(cuni, self.cells_h[can]),
                        sampler.tri_interp_cartuni_lagrange_2d(cuni, self.cells_hu[can]),
                        sampler.tri_interp_cartuni_lagrange_1d(cuni, self.cells_b[can]))
                elif(self.name == 'c'):
                    cbary = sampler.tri_cart_to_bary(np.array(x), self.cells_x[can])
                    # FLASH cell sampling
                    return (sampler.tri_interp_bary_linear_1d(cbary, self.cells_h[can]),
                        sampler.tri_interp_bary_linear_2d(cbary, self.cells_hu[can]),
                        sampler.tri_interp_bary_linear_1d(cbary, self.cells_b[can]))
        return ()

    def __str__(self):
        # Pretty print layer info
        return 'Layer "' + self.name + '" (' + self.layer_name_map[self.name] + \
            '), data layout: ' + str(self.param) + ', HDF5 data ' + \
            ("loaded, " + str(self.cells_length) + " cells" \
            if self.cells_length != 0 else "not yet loaded") + ", bounds: " + str(np.round(self.index_bounds, 4))

class Step:
    """This class represents the metadata of one step in a XMF file, plus methods to manage the data"""

    def __init__(self, index, time, layers):
        self.index = index
        self.time = time
        self.layers = layers

    def num_cells(self):
        csum = 0
        for layer in self.layers:
            csum += layer.cells_length
        return csum

    def load(self):
        for layer in self.layers: layer.load()

    def unload(self):
        for layer in self.layers: layer.unload()

    def sample(self, x):
        for layer in self.layers: 
            res = layer.sample(x)
            if(len(res) != 0): return res
        return (None, (None, None), None)

    def bounds(self):
        xmin = math.inf
        ymin = math.inf
        xmax = -math.inf
        ymax = -math.inf
        for layer in self.layers: 
            xmin = min(xmin, layer.index_bounds[0])
            ymin = min(ymin, layer.index_bounds[1])
            xmax = max(xmax, layer.index_bounds[2])
            ymax = max(ymax, layer.index_bounds[3])
        return (xmin, ymin, xmax, ymax)

    def __str__(self):
        # Pretty print step info
        layerstr = ""
        for i, layer in enumerate(self.layers):
            layerstr += " - [" + str(i) + "] " + str(layer)
            if i < len(self.layers) - 1: layerstr += "\n"
        return 'Computation step #' + str(self.index) + ', t=' + str(self.time) + \
            "s, " + str(len(self.layers)) + " layer(s), " + str(self.num_cells()) + \
            " cells in all layers, bounds: " + str(np.round(self.bounds(), 4)) + "\n" + layerstr

class Reader:
    """This class represents a collection of steps contained in a XMF file, and provides methods to access them"""

    def __init__(self, path):
        self.path = path
        self.steps = []
        self._parse()

    def __str__(self):
        # Pretty print step collection
        return 'XMF file with ' + str(len(self.steps)) + ' steps, loaded from ' + self.path

    def _parse(self):
        tree = etree.parse(self.path)
        root = tree.getroot()

        # Iterate all steps in file
        for grid in root.findall("./Domain/Grid[@GridType='Collection']/Grid"):
            time = float(grid.find("Time").attrib['Value'])
            if time >= 0:
                layers = []

                def subgrid_to_layer(subgrid, lname):
                    # Extract needed data from step XML using Xpath
                    ncells, valstw = list(map(int, subgrid.find("Topology/DataItem").attrib['Dimensions'].split()))
                    h5path, h5g = subgrid.find("Geometry/DataItem").text.split(':')
                    valsgw = int(subgrid.find("Geometry/DataItem").attrib['Dimensions'].split()[1])
                    h5h = subgrid.find("Attribute[@Name='WaterHeight']/DataItem").text.split(':')[1]
                    h5hu = subgrid.find("Attribute[@Name='Momentum']/DataItem").text.split(':')[1]
                    h5b = subgrid.find("Attribute[@Name='Bathymetry']/DataItem").text.split(':')[1]
                    h5l = subgrid.find("Attribute[@Name='Plotter']/DataItem").text.split(':')[1]
                    attrwsym = subgrid.find("Attribute[@Name='WaterHeight']").attrib['Center']
                    # Attr length from XML
                    attrw = 0
                    if attrwsym == "Cell": attrw = 1
                    elif attrwsym == "Node": attrw = valstw
                    # Attr length from HDF5
                    h5path_abs = os.path.join(os.path.dirname(self.path), h5path)
                    h5 = h5py.File(h5path_abs, 'r')
                    attrw_hdf = h5[h5h].shape[1]
                    h5.close()
                    # Add layer object
                    layerparam = Parameter(valsgw, valstw, attrw_hdf, attrw)
                    return StepLayer(layerparam, lname, h5path_abs, h5g, h5h, h5hu, h5b, h5l)

                if ("GridType" in grid.attrib and grid.attrib["GridType"] == "Collection") and \
                        ("CollectionType" in grid.attrib and grid.attrib["CollectionType"] == "Spatial"):
                    # Multi-layer step (ADER-DG)
                    for subgrid in grid.findall("Grid"):
                        layers.append(subgrid_to_layer(subgrid, subgrid.find("Information[@Name='Layer']").attrib['Value']))
                else:
                    # Legacy single-layer (in root grid) step (FLASH)
                    layers.append(subgrid_to_layer(grid, "p"))
                # Add step object
                index = int(grid.find("Attribute[@Name='Step']/DataItem").text)
                self.steps.append(Step(index, time, layers))

    def step_near_time(self, t):
        # Find step closest to timestamp
        if(len(self.steps) != 0):
            times = np.array(list(map(lambda s: s.time, self.steps)))
            return (np.abs(times - t)).argmin()
        else:
            raise Exception("No steps loaded")
