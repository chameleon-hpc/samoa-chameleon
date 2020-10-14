#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import sys, os, getopt

def args(argv):
    inputfile = ''
    outputfile = ''
    scale = 1.0
    sample = 0
    try:
        opts, args = getopt.getopt(argv, 'hi:o:s:a', ['ifile=', 'ofile=', 'scale=', 'sample='])
    except getopt.GetoptError:
        print(sys.argv[0] + ' -i <inputfile> -o <outputfile> -s <scale> -a <sample>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0] + ' -i <inputfile> -o <outputfile> -s <scale> -a <sample>')
            sys.exit(0)
        elif opt in ('-i', '--ifile'):
            inputfile = arg
        elif opt in ('-o', '--ofile'):
            outputfile = arg
        elif opt in ('-s', '--scale'):
            scale = float(arg)
        elif opt in ('-a', '--sample'):
            sample = int(arg)
    return (inputfile, outputfile, scale, sample)

def info(dp):
    yext = (dp.index[0], dp.index[-1])
    ydim = len(dp.index)
    print("Y: " + str(ydim) + " values from " + str(yext[0]) + " to " + str(yext[1]))
    xext = (dp.columns[0], dp.columns[-1])
    xdim = len(dp.columns)
    print("X: " + str(xdim) + " values from " + str(xext[0]) + " to " + str(xext[1]))

def doublesample(dp):
    dfp2 = dp.copy()
    ys_new = [np.nan]*((len(dp.index.values) * 2) - 1)
    for i in range(0, len(dp.index.values)): ys_new[i*2] = dp.index.values[i]
    for y in np.setdiff1d(pd.Series(ys_new).interpolate().values, ys_new): dfp2.loc[y] = np.nan
    xs_new = [np.nan]*((len(dp.columns.values) * 2) - 1)
    for i in range(0, len(dp.columns.values)): xs_new[i*2] = dp.columns.values[i]
    for x in np.setdiff1d(pd.Series(xs_new).interpolate().values, xs_new): dfp2[x] = np.nan
    dfp2 = dfp2.sort_values("y(m)")
    dfp2 = dfp2.sort_values("x(m)", 1)
    dfp2 = dfp2.interpolate('cubic', 0)
    dfp2 = dfp2.interpolate('cubic', 1)
    return dfp2

if __name__ == '__main__':
    inputfile, outputfile, scale, sample = args(sys.argv[1:])
    print("Z Scale: " + str(scale))

    # read csv
    df = pd.read_csv(inputfile)
    dfp = df.pivot(index="y(m)", columns="x(m)", values="z(m)")
    print("Loaded data")
    info(dfp)

    # interpolate two times
    for i in range(0, sample):
        dfp = doublesample(dfp)
        print("Interpolation done "+str(i+1)+" times")
        info(dfp)

    # create netcdf
    root_grp = Dataset(outputfile, 'w', format='NETCDF4')
    root_grp.description = 'Okushiri'

    # dimensions
    root_grp.createDimension('y', len(dfp.index))
    root_grp.createDimension('x', len(dfp.columns))

    # variables
    y = root_grp.createVariable('y', 'f8', ('y',))
    x = root_grp.createVariable('x', 'f8', ('x',))
    z = root_grp.createVariable('z', 'f8', ('y', 'x',))

    # data
    y[:] = dfp.index.values
    x[:] = dfp.columns.values
    z[:,:] = dfp.values * scale
    root_grp.close()
    print("Wrote data")