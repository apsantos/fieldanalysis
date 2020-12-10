#!/usr/bin/python
import math as m
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
    parser.add_argument("--density_collumn", type=int, default=3,
                   help='data column with averaged density [0-indexed]')
    parser.add_argument("--skip_lines", type=int, default=0,
                   help='skip these lines in data file')
    parser.add_argument("--Vfactor", type=float, default=1,
                   help='assume it is the width of the mesh')
    parser.add_argument("--bin_width", type=float, default=0,
                   help='assume it is the width of the mesh')

    if len(parser.parse_args().mesh) != 3:
        print "--mesh must have 3 values for x-, y- and z-direction"
        return
        
    Nx = parser.parse_args().mesh
    L = parser.parse_args().L_direction
    rho_c = parser.parse_args().density_collumn
    start_line = parser.parse_args().skip_lines
    Vfactor = parser.parse_args().Vfactor
    if parser.parse_args().direction == 'x':
        loop_dir = 0
    elif parser.parse_args().direction == 'y':
        loop_dir = 1
    elif parser.parse_args().direction == 'z':
        loop_dir = 2

    if parser.parse_args().bin_width > 0:
        binwidth = parser.parse_args().bin_width
    else:
        binwidth = L/float(Nx[loop_dir])

    Nbin = int(m.ceil(L/binwidth))
    
    # Read in density file
    x, y, z, rho = fts.read_real_all_3d_dat(Nx[0], Nx[1], Nx[2], parser.parse_args().mesh_file, rho_c, start_line)
    
    # Calculate density slice
    if loop_dir == 0:
        slice, std = fts.avg_3d_full_xyz_slice_no_shift(x, loop_dir, Nbin, binwidth, rho)
    elif loop_dir == 1:
        slice, std = fts.avg_3d_full_xyz_slice_no_shift(y, loop_dir, Nbin, binwidth, rho)
    elif loop_dir == 2:
        slice, std = fts.avg_3d_full_xyz_slice_no_shift(z, loop_dir, Nbin, binwidth, rho)
    
    # Write output
    otp = open(parser.parse_args().profile_file,"w")
    
    #dx = x[1] - x[0]
    
    rho0 = parser.parse_args().density
    otp.write("# %s position real_rho real_rho_std imaginary_rho imaginary_rho_std\n" % parser.parse_args().direction)
    for i in range(0, Nbin):
        line = "%lf %e %e %e %e\n" % (i/float(Nbin), Vfactor*slice[i].real/rho0, Vfactor*std[i].real/rho0,
            Vfactor*slice[i].imag/rho0 , Vfactor*std[i].imag/rho0 )
        otp.write(line)
    
    otp.close()

if __name__ == '__main__':
    sys.exit(main())
