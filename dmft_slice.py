#!/usr/bin/python
#import math as m
import fts_utils as fts
import sys, argparse

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Average meshed-field density'
                                                 'and convert to 1-d density profile')
    parser.add_argument("--mesh_file", type=str, required=True,
                   help='input simulation meshed-field density file')
    parser.add_argument("--profile_file", type=str, default='dens_no-shift_slice.dat',
                   help='output averaged density profile file')
    parser.add_argument("--L_direction", type=float, required=True,
                   help='length of simulation box in averaging direction')
    parser.add_argument("--direction", type=str, choices=['x', 'y', 'z'], default='z',
                   help='averaging direction')
    parser.add_argument("--mesh", type=int, nargs='+', required=True,
                   help='mesh size in the x-,y- and z-directions.')
    parser.add_argument("--density", type=float, required=True,
                   help='density, rho_0, in field part of simulation')

    if len(parser.parse_args().mesh) != 3:
        print "--mesh must have 3 values for x-, y- and z-direction"
        return
        
    Nx = parser.parse_args().mesh
    L = parser.parse_args().L_direction
    if parser.parse_args().direction == 'x':
        loop_dir = 0
    elif parser.parse_args().direction == 'y':
        loop_dir = 1
    elif parser.parse_args().direction == 'z':
        loop_dir = 2
    
    # Read in density file
    x, y, z, rho = fts.read_real_3d_dat(Nx[0], Nx[1], Nx[2], parser.parse_args().mesh_file)
    
    # Calculate density slice
    slice, std = fts.avg_3d_full_slice_no_shift(Nx[0], Nx[1], Nx[2], loop_dir, rho)
    
    # Write output
    otp = open(parser.parse_args().profile_file,"w")
    
    #dx = x[1] - x[0]
    dx = L / Nx[loop_dir]
    
    rho0 = parser.parse_args().density
    otp.write("# %s position real_rho real_rho_std imaginary_rho imaginary_rho_std\n" % parser.parse_args().direction)
    for i in range(0, Nx[loop_dir]):
        line = "%lf %e %e %e %e\n" % (i*dx/L, slice[i].real/rho0, std[i].real/rho0 , 
            slice[i].imag/rho0 , std[i].imag/rho0 )
        otp.write(line)
    
    otp.close()

if __name__ == '__main__':
    sys.exit(main())
