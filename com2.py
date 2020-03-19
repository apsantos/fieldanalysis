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

    ndat = file_len(parser.parse_args().profile_file) - parser.parse_args().start_line
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
        density[i,0] = float(data[0])
        density[i,1] = float(data[1])
        i += 1
 
    den_sum = np.sum(density[:,1])
  
    error = 100000
    tol = 1e-6
    while (abs(error) > tol):
        # calculate the center-of-mass of the density profile
        COM = 0.0
        for i in range(ndat):
            COM += density[i,1] * (density[i,0] - m.floor(density[i,0]))
        COM = COM / den_sum
        error = COM - 0.5
      
        # shift the profile by the center-of-mass
        for i in range(ndat):
            density[i,0] = density[i,0] - error - m.floor(density[i,0] - error)

    # sort the values because they were shifted
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
