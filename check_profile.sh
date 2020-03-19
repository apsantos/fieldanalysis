#!/bin/bash
# Process the simulation density profile fits
# Assess if it is big and equilibrated enough to be used in a phase diagram

read_input=1
input_file=$1
if [ $read_input == 1 ]; then
    sigma=`sed '7!d' $input_file | awk '{print $1}'`
    NG=`sed '7!d' $input_file | awk '{print $2}'`
    Rp=`sed '8!d' $input_file | awk '{print $1}'`
    rho=`sed '13!d' $input_file | awk '{print $1}'`
    Lx=`sed '18!d' $input_file | awk '{print $1}'`
    Ly=`sed '18!d' $input_file | awk '{print $2}'`
    Lz=`sed '18!d' $input_file | awk '{print $3}'`
    Mx=`sed '19!d' $input_file | awk '{print $1}'`
    My=`sed '19!d' $input_file | awk '{print $2}'`
    Mz=`sed '19!d' $input_file | awk '{print $3}'`
    phiP=`sed '6!d' $input_file | awk '{print $1}'`
    Vp=`python -c 'print('${rho}'*12.56637*'${Rp}'**3.0/3.0)'`
    NgNP=`python -c 'print(int('$phiP'*'$Lx'*'$Ly'*'$Lz'*'$rho'/'$Vp'))'`
else
    sigma=$2
    Rp=$3
    NG=$4
    rho=$5
    Mx=$6
    My=$7
    Mz=$8
    Lx=$9
    Ly=${11}
    Lz=${12}
fi

# get profile fit data
avg_rhog=`tail -n1 avg_rhodb.fits | awk '{print $1}'`
avg_rhol=`tail -n1 avg_rhodb.fits | awk '{print $2}'`
avg_rho=`python -c 'print(max('$avg_rhog','$avg_rhol'))'`
avg_xi=`tail -n1 avg_rhodb.fits | awk '{print $3}'`
avg_z0=`tail -n1 avg_rhodb.fits | awk '{print $4}'`
rhog=`tail -n1 rhodb.fits | awk '{print $1}'`
rhol=`tail -n1 rhodb.fits | awk '{print $2}'`
rho=`python -c 'print(max('$rhog','$rhol'))'`
xi=`tail -n1 rhodb.fits | awk '{print $3}'`
z0=`tail -n1 rhodb.fits | awk '{print $4}'`

# is the phase well-developed?
fraction1=`python -c 'print('$avg_z0'-3*'$avg_xi')'`
fraction2=`python -c 'print('$avg_z0'+3*'$avg_xi')'`

# are there enough NP?
dilutemin=`min_gNP.py --direction z --L $Lx $Ly $Lz --fitted_phiM $avg_rho --fitted_width $avg_xi --fitted_z0 $avg_z0 --NgNP $NgNP --Rp $Rp --sigma $sigma --NG $NG | awk '{print $3}'`
simdilute=`python -c 'print(1.0-'$avg_rho')'`

# is it equilibrated?
Dphi_gNP=`python -c 'print(('$avg_rhog'-'$rhog')/'$avg_rhog')'`
Dphi_M=`python -c 'print(('$avg_rho'-'$rho')/'$avg_rho')'`
Dxi=`python -c 'print(('$avg_xi'-'$xi')/'$avg_xi')'`
Dz0=`python -c 'print(('$avg_z0'-'$z0')/'$avg_z0')'`

#echo "# developed?, enoughNP?, equilibrated?"
#echo "# (z0-3Delta) (z0+3Delta) phigNP-low minphi D(phi_gNP) D(phi_M) D(delta) D(z0)" 
echo $fraction1 $fraction2 $dilutemin $simdilute $Dphi_gNP $Dphi_M $Dxi $Dz0
