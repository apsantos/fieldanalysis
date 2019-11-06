#for i in 27; do for j in 0  0.25  0.5  0.75  1   ; do cd ${i}/${j}/; dmft_slice.py 25 25 315 250 2 avg_rhodb.dat 3.94342; com.py dens_no-shift_slice.dat tmp.dat; fit_curve.py 1 tmp.dat tmp2.dat 1; cd ..; cd ..; done; done
for i in 27; do 
    for j in 0  0.25  0.5  0.75  1; do 
        cd ${i}/${j}
        dmft_slice.py 25 25 315 250 2 avg_rhodb.dat 3.94342
        com.py dens_no-shift_slice.dat tmp.dat 
        fit_curve.py 1 tmp.dat tmp2.dat 1
        cd ../..
    done
done
