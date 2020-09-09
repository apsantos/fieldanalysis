#!/usr/bin/python
import sys, argparse
import math as m
import os, sys
import numpy as np
"""
Algorithm from Jason
"""

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Calculate radial distriution from density files')
    parser.add_argument('-i', "--grid_file", type=str, required=True,
                   help='grid data file. with the following format:'
                        'ix iy iz rho_1 rho_2 .. rho_3. where rho_i is density of a component')
    parser.add_argument('-o', "--density_file", type=str, help='output file name')
    parser.add_argument("--start_line", type=int, default=1,
                   help='Set the starting read in line. Default, 1 (skip first line).')
    parser.add_argument('-r', "--np_radius", type=float, required=True, help='np radius')
    parser.add_argument('-n', "--c_npdensity", type=int, default=4,
                   help='column with the NP/center of radial density to calculate profile')
    parser.add_argument('-d', "--c_density", type=int, default=3,
                   help='column with the density to calculate profile')
    parser.add_argument('-b', "--density_binwidth", type=float, required=True,
                   help='Bin width')

    ifile = open(parser.parse_args().grid_file, 'r')
    if parser.parse_args().density_file:
        ofile = open(parser.parse_args().density_file, 'w')
    else:
        ofile = open('ouput.rho', 'w')
    start_line = parser.parse_args().start_line
    c_density = parser.parse_args().c_density
    c_npdensity = parser.parse_args().c_npdensity
    binwidth = parser.parse_args().density_binwidth
    Rnp = parser.parse_args().np_radius

    ngrid = file_len(parser.parse_args().grid_file)
    ngrid -= start_line

    pos = np.zeros( (ngrid,3) )
    rho = np.zeros( (ngrid) )
    rhonp = np.zeros( (ngrid) )

    for i in range(start_line):
        ifile.readline().split()

    for i in range(ngrid):
        line = ifile.readline().split()
        pos[i,0] = float( line[0] )
        pos[i,1] = float( line[1] )
        pos[i,2] = float( line[2] )
        rho[i] = float( line[c_density] )
        rhonp[i] = float( line[c_npdensity] )
      
    ifile.close()


    Lmax = np.zeros( 3 )
    Lmax[0] = np.max( pos[:,0] )
    Lmax[1] = np.max( pos[:,1] )
    Lmax[2] = np.max( pos[:,2] )
    Lmin = np.zeros( 3 )
    Lmin[0] = np.min( pos[:,0] )
    Lmin[1] = np.min( pos[:,1] )
    Lmin[2] = np.min( pos[:,2] )

    nbins = int( np.max( Lmax ) / binwidth ) 
    
    norm = np.zeros( nbins, 'd' )
    rad = np.zeros( nbins , 'd' )

    # find the NP center(s)
    cent = (Lmax - Lmin) / 2.0 

    maxrhonp = max(rhonp)
    #find centers
    center_grid = []
    for i in range(ngrid):
        if rhonp[i] == maxrhonp:
            center_grid.append(i)

    maxpos_south = Lmin
    minpos_north = Lmax
    maxpos = Lmin
    minpos = Lmax
    if len(center_grid) == 1:
        cent = pos[center_grid[0],:]
    else:
        for ig in range(len(center_grid)):
            for idim in range(3):
                maxpos[idim] = max(maxpos[idim], pos[ig, idim])
                minpos[idim] = min(minpos[idim], pos[ig, idim])
                if pos[ig, idim] < cent[idim]:
                    maxpos_south[idim] = max(maxpos_south[idim], pos[ig, idim])
                else:
                    minpos_north[idim] = min(minpos_north[idim], pos[ig, idim])

        print "particle position: " , cent
        print maxpos, minpos, maxpos_south, minpos_north
        maxd = 2.5*Rnp
        for idim in range(3):
            # split by boundary
            if maxpos[idim] == Lmax[idim] and minpos[idim] == Lmin[idim]:
                adiff = Lmax[idim] - minpos_north[idim]
                bdiff = maxpos_south[idim] - Lmin[idim]
                if adiff > bdiff:
                    cent[idim] = minpos_north[idim] + bdiff
                else:
                    cent[idim] = maxpos_south[idim] + adiff
                
            else:
                cent[idim] = (maxpos[idim] - minpos[idim])/2.0
        
    print "particle position: " , cent

    # radial density profile
    for i in range(ngrid):
        mdr2 = 0.0
        for j in range(3):
            dr = pos[i,j] - cent[j]
            mdr2 += (pos[i,j] - cent[j]) ** 2

        mdr = m.sqrt( mdr2 )

        cur_bin = int( mdr / binwidth )

        if ( cur_bin < nbins ):
            norm[ cur_bin ] += 1
            rad[ cur_bin ] += rho[i]

    ofile.write( "# r rho(r)" )
    
    for i in range(nbins):
        if ( norm[i] > 0.0 ):
            line = '%f %f\n' % ( (i+0.5) * binwidth , rad[i]/norm[i] )
        else:
            line = '%f 0.0\n' % ( (i+0.5) * binwidth )

        ofile.write( line )

    ofile.close()

if __name__ == '__main__':
    sys.exit(main())
