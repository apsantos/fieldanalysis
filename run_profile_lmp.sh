#!/bin/bash
# Run scripts to get information on the density profiles
read_input=$1
echo $input_file
if [ $read_input == 1 ]; then
    input_file=$2
    rho=`grep "actual rho" $input_file | tail -n1 | awk '{print $9}'`
    Mx=`grep "mesh" $input_file | tail -n1 | awk '{print $3}'`
    My=`grep "mesh" $input_file | tail -n1 | awk '{print $4}'`
    Mz=`grep "mesh" $input_file | tail -n1 | awk '{print $5}'`

    datafile=`grep "read_data" $input_file | tail -n1 | awk '{print $2}'`
    Lx=`head $datafile | grep xlo | awk '{print $2}'`
    Ly=`head $datafile | grep ylo | awk '{print $2}'`
    Lz=`head $datafile | grep zlo | awk '{print $2}'`

    Rp=`grep "tild/coeff" $input_file | grep "erfc" | tail -n1 | awk '{print $6}'`
    Vfactor=`python -c 'print(4./3.*3.14159*'$rho'*'$Rp'*'$Rp'*'$Rp')'`
    avgfile=`grep "tild/ave/grid" $input_file | tail -n1 | awk '{print $6}'`
    instfile=`grep "tild/write_grid_data" $input_file | tail -n1 | awk '{print $4}'`

else
    rho=$2
    Mx=$3
    My=$4
    Mz=$5
    Lz=$6
    instfile=$7

fi


echo $instfile $Mx $My $Mz $Lz $rho 
# average the density mesh and create a profile
dmft_slice.py --mesh_file $instfile --profile_file rho_1.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho --density_collumn 3 --skip_lines 1
#dmft_slice.py --mesh_file $avgfile --profile_file avgrho_1.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho --density_collumn 3 --skip_lines 1
dmft_slice.py --mesh_file $instfile --profile_file rho_2.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho --density_collumn 4 --skip_lines 1
#dmft_slice.py --mesh_file $avgfile --profile_file avgrho_2.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho --density_collumn 4 --skip_lines 1
dmft_slice.py --mesh_file $instfile --profile_file rho_3.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho --density_collumn 5 --skip_lines 1 --Vfactor $rho #--bin_width 4
#dmft_slice.py --mesh_file $instfile --profile_file rho_3.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho --density_collumn 5 --skip_lines 1 --Vfactor $Vfactor #--bin_width 4
#dmft_slice.py --mesh_file $avgfile --profile_file avgrho_3.profile --L_direction $Lz --direction z --mesh $Mx $My $Mz --density $rho --density_collumn 5 --skip_lines 1

# center the denisty profile
com2.py --profile_file rho_1.profile --centered_file rho_1.shifted -l 1
#com2.py --profile_file avg_rho_1.profile --centered_file avg_rho_1.shifted -l 1
com2.py --profile_file rho_2.profile --centered_file rho_2.shifted -l 1
#com2.py --profile_file avg_rho_2.profile --centered_file avg_rho_2.shifted -l 1
#com2.py --profile_file rho_3.profile --centered_file rho_3.shifted -l 1
#com2.py --profile_file avg_rho_3.profile --centered_file avg_rho_3.shifted -l 1

# fit the denisty profile
fit_curve.py --profile_file  rho_1.shifted --fit_file rho_1.shifted_fit --L_direction 1 -l 1 --fit_method odr > rho_1.fits
#fit_curve.py --profile_file  avg_rho_1.shifted --fit_file avg_rho_1.shifted_fit --L_direction 1 -l 1 --fit_method odr > avg_rho_1.fits 
fit_curve.py --profile_file  rho_2.shifted --fit_file rho_2.shifted_fit --L_direction 1 -l 1 --fit_method odr > rho_2.fits
#fit_curve.py --profile_file  avg_rho_2.shifted --fit_file avg_rho_2.shifted_fit --L_direction 1 -l 1 --fit_method odr > avg_rho_2.fits 
#fit_curve.py --profile_file  rho_3.shifted --fit_file rho_3.shifted_fit --L_direction 1 -l 1 --fit_method odr > rho_3.fits
#fit_curve.py --profile_file  avg_rho_3.shifted --fit_file avg_rho_3.shifted_fit --L_direction 1 -l 1 --fit_method odr > avg_rho_3.fits 

