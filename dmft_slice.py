#!/usr/bin/python
import numpy as nmp
import os, sys
import cmath as cm
import math as m
import fts_utils as fts

if (len(sys.argv) < 7):
  print "Usage: 3D_dens_slice.py [Nx] [Ny] [Nz] [L-direction] [direction] [density file] [rho0]"
  exit(1)

Nx = nmp.zeros(3,'int')
Nx[0] = int( sys.argv[1] )
Nx[1] = int( sys.argv[2] )
Nx[2] = int( sys.argv[3] )
L = float(sys.argv[4])
dirs = nmp.zeros(3, 'int')
loop_dir = int( sys.argv[5] )
#notloop = nmp.where(dirs != loop_dir)[0]
#Mnl = Nx[ notloop[0] ] * Nx[ notloop[1] ]

# Read in density file
x, y, z, rho = fts.read_real_3d_dat(Nx[0], Nx[1], Nx[2], sys.argv[6])

# Calculate density slice
slice , std = fts.avg_3d_full_slice_no_shift(Nx[0], Nx[1], Nx[2], loop_dir, rho)

# Write output
otp = open("dens_no-shift_slice.dat","w")

#dx = x[1] - x[0]
dx = L / Nx[loop_dir]

rho0 = float(sys.argv[7])
for i in range(0, Nx[loop_dir]):
  line = "%lf %e %e %e %e\n" % (i*dx/L, slice[i].real/rho0, std[i].real/rho0 , 
      slice[i].imag/rho0 , std[i].imag/rho0 )
  otp.write(line)

otp.close()
