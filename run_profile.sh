#!/bin/bash
# Run scripts to get information on the density profiles
read_input=1
input_file=$1
if [ $read_input == 1 ]; then
    rho=`sed '13!d' $input_file | awk '{print $1}'`
    Lx=`sed '18!d' $input_file | awk '{print $1}'`
    Ly=`sed '18!d' $input_file | awk '{print $2}'`
    Lz=`sed '18!d' $input_file | awk '{print $3}'`
    Mx=`sed '19!d' $input_file | awk '{print $1}'`
    My=`sed '19!d' $input_file | awk '{print $2}'`
    Mz=`sed '19!d' $input_file | awk '{print $3}'`
else
    rho=$5
    Mx=$6
    My=$7
    Mz=$8
    Lx=$9
    Ly=${11}
    Lz=${12}
fi

# average the density mesh and create a profile
dmft_slice.py --mesh_file rhodb.dat --profile_file rhodb.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho
dmft_slice.py --mesh_file avg_rhodb.dat --profile_file avg_rhodb.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho

# center the denisty profile
com2.py --profile_file rhodb.profile --centered_file rhodb.shifted -l 1
com2.py --profile_file avg_rhodb.profile --centered_file avg_rhodb.shifted -l 1

# fit the denisty profile
fit_curve.py --profile_file  rhodb.shifted --fit_file rhodb.shifted_fit --L_direction 1 -l 1 --fit_method odr > rhodb.fits 
fit_curve.py --profile_file  avg_rhodb.shifted --fit_file avg_rhodb.shifted_fit --L_direction 1 -l 1 --fit_method odr > avg_rhodb.fits 

