#!/usr/bin/python
import sys, argparse
import numpy as np
from scipy import integrate

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit vapor, liquid density to get critical point and curve')
    parser.add_argument("--profile_file", type=str, required=True, help='input density profile file')
    parser.add_argument('-l', "--start_line", type=int, default=1, help='starting line of density profile file')
    parser.add_argument("--r_col", type=int, default=0, help='radius data collumn')
    parser.add_argument("--r_min", type=float, default=0.0, help='radius min to integrate')
    parser.add_argument("--r_max", type=float, default=1000000.0, help='radius min to integrate')
    parser.add_argument("--rho_col", type=int, default=1, help='density data collumn')
    parser.add_argument("--rho_inf", type=float, default=1.0, help='the far-field density')

    c_r = parser.parse_args().r_col
    c_rho = parser.parse_args().rho_col
    rho_inf = parser.parse_args().rho_inf
    r_min = parser.parse_args().r_min
    r_max = parser.parse_args().r_max
    ndat = file_len(parser.parse_args().profile_file) - parser.parse_args().start_line

    # read in temeprature and vapor/liquid densities
    inp = open(parser.parse_args().profile_file, 'r') 

    j = 0
    r = []
    deltarho = []
    for line in inp:
        data = line.strip().split()
        if j < parser.parse_args().start_line:
            j += 1
            continue
        tr = float(data[c_r])
        if r_min < tr < r_max:
            r.append( tr )
            deltarho.append( float(data[c_rho])-rho_inf )

    r = np.array(r)
    deltarho = np.array(deltarho)
    integral = integrate.simps(deltarho,x=r)
    
    print integral
if __name__ == '__main__':
    sys.exit(main())

