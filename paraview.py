#!/usr/bin/python

import sys, argparse

    
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Write paraview file')
    parser.add_argument("--grid_file", type=str, required=True, help='input density profile file')
    parser.add_argument("--tec_file", type=str, default='rho.tec', help='output fitted density profile file')
    parser.add_argument("--Lx", type=float, required=True, help='length of simulation box in x direction')
    parser.add_argument("--Ly", type=float, required=True, help='length of simulation box in y direction')
    parser.add_argument("--Lz", type=float, required=True, help='length of simulation box in z direction')
    parser.add_argument('-l', "--start_line", type=int, default=1, help='starting line of density profile file')
    parser.add_argument('-c', "--rho_col", type=int, default=3, help='column with density')

    Lx = parser.parse_args().Lx
    Ly = parser.parse_args().Ly
    Lz = parser.parse_args().Lz
    col = parser.parse_args().rho_col

    out = open(parser.parse_args().tec_file, 'w')
    out.write("TITLE = \"%s\"\n" % parser.parse_args().tec_file)
    out.write("VARIABLES = \"X\", \"Y\", \"Z\", \"Real\"" ); 
    out.write("ZONE I=%d, J=%d, K=%d, F=POINT\n" % (Lx, Ly, Lz) ); 
#
    # read in position and density
    inp = open(parser.parse_args().grid_file, 'r') 

    j = 0
    for line in inp:
        if j < parser.parse_args().start_line:
            j += 1
            continue
        data = line.strip().split()
        x = float(data[0])
        y = float(data[1]) 
        z = float(data[2]) 
        rho = float(data[col]) 
        out.write("%f %f %f %f\n" % (x, y, z, rho) ); 

    out.close()
    inp.close()
if __name__ == '__main__':
    sys.exit(main())
