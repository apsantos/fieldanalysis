#!/bin/bash
rho=$1
dmft_slice.py 25 25 315 250 2 avg_rhodb.dat $rho
com2.py dens_no-shift_slice.dat tmp.dat 
fit_curve.py 1 tmp.dat tmp2.dat 1
