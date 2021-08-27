#!/usr/bin/python
import sys, argparse
import math as m
import numpy as np

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='shift 1-d density profile so that the dense part is in the middle')
    parser.add_argument("--profile_file", type=str, required=True, help='input density profile file')
    parser.add_argument("--centered_file", type=str, default='centered_profile.dat', help='output centered density profile file')
    parser.add_argument('-l', "--start_line", type=int, default=1, help='starting line of density profile file')
    parser.add_argument("--r_collumn", type=int, default=0, help='collumn for distance')
    parser.add_argument("--rho_collumn", type=int, default=1, help='collumn for density')
    parser.add_argument("--rho_normalizer", type=float, default=1.0, help='normalize density to [0,1]')
    parser.add_argument("--r_normalizer", type=float, default=1.0, help='normalize distance to [0,1]')

    ndat = file_len(parser.parse_args().profile_file) - parser.parse_args().start_line
    r_c = parser.parse_args().r_collumn
    rho_c = parser.parse_args().rho_collumn
    r_n = parser.parse_args().r_normalizer
    rho_n = parser.parse_args().rho_normalizer
    density = np.zeros((ndat,2),'d')

    # read in box position and density
    inp = open(parser.parse_args().profile_file, 'r') 
    i = 0
    j = 0
    for line in inp:
        data = line.strip().split()
        if j < parser.parse_args().start_line:
            j += 1
            continue
        density[i,0] = float(data[r_c]) / r_n
        density[i,1] = float(data[rho_c]) / rho_n
        i += 1
 
    den_sum = np.sum(density[:,1])
  
    error = 100000
    tol = 1e-7
    firsttime = 1
    totalerror = 0
    while (abs(error) > tol):
        # calculate the center-of-mass of the density profile
        COM = 0.0
        for i in range(ndat):
            COM += density[i,1] * (density[i,0] - m.floor(density[i,0]))
        COM = COM / den_sum
        error = COM - 0.5
        totalerror += error
      
        # shift the profile by the center-of-mass
        for i in range(ndat):
            density[i,0] = density[i,0] - error - m.floor(density[i,0] - error)

    # sort the values because they were shifted
    print "COM/Lz: ", totalerror+0.5
    density = density[np.argsort(density[:,0])] 

    # center the dense part
    if density[0,1] > density[int(ndat/2.0),1]:
        for i in range(ndat):
            if density[i,0] >= 0.5: 
                density[i,0] -= 0.5
            elif density[i,0] < 0.5: 
                density[i,0] += 0.5

    # sort the values because they were shifted
    density = density[np.argsort(density[:,0])] 
  
    otp = open(parser.parse_args().centered_file, 'w') 
    otp.write("# position density\n")
    for i in range(ndat):
        line = '%lf %lf\n' % (density[i,0], density[i,1])
        otp.write(line)
  
    otp.close()

if __name__ == '__main__':
    sys.exit(main())
