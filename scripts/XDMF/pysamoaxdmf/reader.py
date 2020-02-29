# Samoa^2 XDMF Reader for Python 3

# This module provides a easy-to-use object-orientated interface
# for the XDMF (HDF5 database and XML-based XMF index) data generated
# by the XDMF output module in samoa^2

import numpy as np
import h5py
import os
from lxml import etree

class Parameter:
    """This class represents the XDMF data and file layout configuration"""

    def __init__(self, valsg_width, valst_width, attr_width):
        self.valsg_width = valsg_width
        self.valst_width = valst_width
        self.attr_width = attr_width

class Step:
    """This class represents the metadata of one step in a XMF file, plus a HDF5 file handle to the data, and methods to access this data"""

    def __init__(self, index, time, \
            h5path, h5gname, h5hname, h5huname, h5bname):
        self.index = index
        self.time = time
        self.h5path = h5path
        self.h5gname = h5gname
        self.h5hname = h5hname
        self.h5huname = h5huname
        self.h5bname = h5bname
        self.cells_length = 0
        self.unload()

    def __str__(self):
        # Pretty print step info
        return 'Computation step #' + str(self.index) + ', t=' + str(self.time) + \
            "s, HDF5 data " + ("loaded, " + str(self.cells_length) + " cells" if self.cells_length != 0 else "not yet loaded")

    def load(self):
        h5 = h5py.File(self.h5path, 'r')
        h5dg = h5[self.h5gname]
        h5dh = h5[self.h5hname]
        h5dhu = h5[self.h5huname]
        h5db = h5[self.h5bname]
        # Actually load the data to speed up lookups
        gdata = h5dg[...]
        self.cells_length = h5dh.shape[0]
        self.cells_h = h5dh[...]
        self.cells_hu = h5dhu[...]
        self.cells_b = h5db[...]
        # Iterate all cells and unpack location
        self.cells_x = np.empty([self.cells_length, h5dh.shape[1], gdata.shape[1]]) 
        for i in range(0, self.cells_length):
            self.cells_x[i] = gdata[(i * h5dh.shape[1]):(i * h5dh.shape[1]) + h5dh.shape[1]]
        h5.close()   

    def unload(self):
        self.cells_length = 0
        self.cells_x = np.array([])
        self.cells_h = np.array([])
        self.cells_hu = np.array([])
        self.cells_b = np.array([])

class Reader:
    """This class represents a collection of steps contained in a XMF file, and provides methods to access them"""

    def __init__(self, path, param, subgroup):
        self.path = path
        self.param = param
        self.subgroup = subgroup
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
                subgroup_found = False
                for subgrid in root.findall("Grid"):
                    if subgrid.find("Information[@Name='Layer']").attrib['Value'] == self.subgroup:
                        subgroup_found = True
                        # Extract needed data from step XML using Xpath
                        index = int(subgrid.find("Attribute[@Name='Step']/DataItem").text)
                        ncells, valstw = list(map(int, subgrid.find("Topology/DataItem").attrib['Dimensions'].split()))
                        h5path, h5g = subgrid.find("Geometry/DataItem").text.split(':')
                        valsgw = int(subgrid.find("Geometry/DataItem").attrib['Dimensions'].split()[1])
                        h5h = subgrid.find("Attribute[@Name='WaterHeight']/DataItem").text.split(':')[1]
                        h5hu = subgrid.find("Attribute[@Name='Momentum']/DataItem").text.split(':')[1]
                        h5b = subgrid.find("Attribute[@Name='Bathymetry']/DataItem").text.split(':')[1]
                        attrwsym = subgrid.find("Attribute[@Name='WaterHeight']").attrib['Center']
                        attrw = 0
                        if attrwsym == "Cell": attrw = 1
                        elif attrwsym == "Node": attrw = valstw

                        # Sanity check layout
                        if self.param.valsg_width != valsgw:
                            raise Exception("Layout parameter mismatch, geometry should have " + str(self.param.valsg_width) + " spatial dimensions, but has " + str(valsgw))
                        if self.param.valst_width != valstw:
                            raise Exception("Layout parameter mismatch, cell should have " + str(self.param.valst_width) + " geometry entries, but has " + str(valstw))
                        if self.param.attr_width != attrw:
                            raise Exception("Layout parameter mismatch, cell should have " + str(self.param.attr_width) + " attribute entries, but has " + str(attrw))

                        # Add step object
                        self.steps.append(Step(index, time, os.path.join(os.path.dirname(self.path), h5path), h5g, h5h, h5hu, h5b))
                if not subgroup_found:
                    raise Exception("Could not find subgroup '" + self.subgroup + "')

    def step_near_time(self, t):
        # Find step closest to timestamp
        if(len(self.steps) != 0):
            times = np.array(list(map(lambda s: s.time, self.steps)))
            return (np.abs(times - t)).argmin()
        else:
            raise Exception("No steps loaded")
