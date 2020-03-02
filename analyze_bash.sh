#!/bin/bash
rho=$1
Mx=$2
My=$3
Mz=$4
Lz=$5
dmft_slice.py $Mx $My $Mz $Lz 2 avg_rhodb.dat $rho
com2.py dens_no-shift_slice.dat tmp.dat 
fit_curve.py 1 tmp.dat tmp2.dat 1 > fits 
