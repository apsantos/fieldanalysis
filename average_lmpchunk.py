#!/usr/bin/python
import sys, argparse
import math as m
import os, sys
import numpy as np

def getfileinfo(fname):
    f = open(fname, 'r')
    f.readline()
    f.readline()
    f.readline()
    nsteps = 1
    line = f.readline()
    if len(line.strip().split()) > 2:
        return -1, -1
    nchunk = int(line.strip().split()[1])
    ichunk = 0
    maxnchunk = 0
    while line:
        line = f.readline()
        if ichunk < nchunk:
            ichunk += 1
        elif ichunk == nchunk:
            if len(line.strip().split()) < 2: break
            ichunk = 0
            nsteps += 1
            nchunk = int(line.strip().split()[1])
            maxnchunk = max(maxnchunk, nchunk)
        
    maxnchunk = max(maxnchunk, nchunk)
    f.close()

    return nsteps, maxnchunk

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Calculate histogram, average, uncertainty and distirbution from lammps')
    parser.add_argument('-i', "--chunk_file", type=str, required=True,
                   help='LAMMPS ave/time output file from chunk.')
    parser.add_argument('-a', "--average", action="store_true", default=True, help='calculate averages')
    parser.add_argument('-d', "--std_dev", action="store_true", default=True, help='calculate standrd deviations')
    parser.add_argument('-b', "--histogram", type=int, default=0, help='histogram number of bins, set to 0 to ignore')
    parser.add_argument("--calctrace", action="store_true", default=True, help='calctrace')

    ifile = open(parser.parse_args().chunk_file, 'r')
    nsteps, maxnchunk = getfileinfo(parser.parse_args().chunk_file)
    if nsteps == -1:
        print parser.parse_args().chunk_file, " is not a LAMMPS chunk file, ending"
        return    

    ifile.readline(); ifile.readline()
    nbins = parser.parse_args().histogram
    cols = ifile.readline().strip().split()[2:]
    ncol = len(cols)
    istep = 0
    calctrace = parser.parse_args().calctrace
    if calctrace:
        tdata = np.zeros((nsteps,ncol+1))
        tncol = ncol + 1
    else:
        tdata = np.zeros((nsteps,ncol))
        tncol = ncol 
 
    if parser.parse_args().average:
        afile = open(parser.parse_args().chunk_file+'avg', 'w')
        afile.write('# step avg,std ')
        for ic in range(ncol):
            afile.write('%s ' % cols[ic])
        if calctrace:
            afile.write('trace')
        afile.write('\n')

    skipstring = ""
    for ic in range(ncol):
        skipstring += " 0"

    for istep in range(nsteps):
        nchunk = int(ifile.readline().strip().split()[1])
        if nchunk <= 1: continue
        data = []
        for ic in range(tncol):
            data.append( [] )

        for irow in range(nchunk):
            rawdata = ifile.readline()
            if skipstring not in rawdata:
                strdata = rawdata.strip().split()
                for ic in range(ncol):
                    data[ic].append( float(strdata[ic+1]) )
                if calctrace:
                    data[ncol].append( float(strdata[1]) + float(strdata[2]) + float(strdata[3]) )

        data = np.array(data).T

        if parser.parse_args().average:
            afile.write('%d ' % istep)
            for ic in range(tncol):
                afile.write('%f %f ' % (np.mean(data[:,ic]), np.std(data[:,ic])) )
            afile.write('\n')
        if nbins > 0:
            for ic in range(tncol):
                hfile = open(parser.parse_args().chunk_file+'histo'+str(ic), 'w')
                hfile.write('# bin value rho\n')
                histo,edges = np.histogram(data[:,ic],bins=nbins, density=True)
                for i in range(nbins):
                    hfile.write( '%f %f\n' % (edges[i], histo[i]) )
                hfile.close()

    afile.close

if __name__ == '__main__':
    sys.exit(main())
