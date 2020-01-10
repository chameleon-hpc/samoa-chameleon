# This module computes the L2 and Lsup limiter errors given an analytical solution

import sys, os, getopt, psutil

def args(argv):
    inputfile = ''
    datadir = ''
    outputfile = ''
    scenario = ''
    workers = 1
    try:
        opts, args = getopt.getopt(argv, 'hi:d:o:s:n:', ['ifile=', 'ddir=', 'ofile=', 'scenario=', 'nworkers='])
    except getopt.GetoptError:
        print(sys.argv[0] + ' -i <inputfile> -o <outputfile> -d <datadir> -s <scenario> -n <workers>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0] + ' -i <inputfile> -o <outputfile> -d <datadir> -s <scenario> -n <workers>')
            sys.exit(0)
        elif opt in ('-i', '--ifile'):
            inputfile = arg
        elif opt in ('-d', '--ddir'):
            datadir = arg
        elif opt in ('-o', '--ofile'):
            outputfile = arg
        elif opt in ('-s', '--scenario'):
            scenario = arg
        elif opt in ('-n', '--nworkers'):
            workers = int(arg)
    return (inputfile, datadir, outputfile, scenario, workers)

def siginth(sig, frame):
    # Ensure all workers are terminated on interrupt
    print("Process #" + str(os.getpid()) + " got interrupt, terminating workers")
    parent = psutil.Process(os.getpid())
    for child in parent.children(): 
        try:
            child.kill()
        except:
            pass
    sys.exit(1)