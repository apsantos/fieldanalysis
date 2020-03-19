# Scripts to process and analyze TILD simulation output

## Calculate the density profiles in the z direction and fit an error function to the two interfaces:
### run_profile.sh:
Bash script that runs the key python scripts
### check_profile.sh:
Bash script that checks the quality of density profile phase-coexistence predictions; averages and gets uncertainties
### min_gNP.py:
Python script to determine the minimum volume fraction of grafted nanoparticles in the dilute region based on box size
-- `min_gNP.py -h` for info
### dmft_slice.py: 
Python script that takes an output density file and averages over a chosen plane.
-- `dmft_slice.py -h` for info
-- Nx = Number of grid points in x-direction
-- Lx = Length of box in x-direction
-- Hardcoded output is always 'dens_no-shift_slice.dat'
-- Modified from Jason Koski's script.

### com2.py: 
Shifts input file so that the density file is centered at Lz/2

- com.py [input file, typically 'dens_no-shift_slice.dat'] [centered output density]

- com2.py is basically the same file as com.py but it shifts the density file to be centered at 0
-- `com2.py -h` for info
-- Modified from Jason Koski's script.

Once in a while in doing my mass analysis, some of the density profiles would be shifted so that the matrix-dense region was more dense in the center of the box. 'fit_curve.py' is best setup to calculate density profiles for [matrix-dense, gnp-dense, matrix-dense] curves. This way, you have both centering scripts incase you run into this problem.

### fit_curve.py: fits a density curve to 2 error functions and spits out the matrix-dense volume fraction and the gnp-dense volume fraction
--Note dmft_slice.py spits out the density file as phi vs z/Lz, so always use Lz is 1
-- `fit_curve.py -h` for info
--If you get funky results from this, you likely have to center the density curve using com2.py OR once in a blue moon, you may have to change your initial guess.
--Note you’ll need the scripts in the fts folder as well, as they have some functions that are used.
-- Modified from Jason Koski's script.

## Trajectory analysis
### np_properties.py
Analyzes the trajectories. 
-- `np_properties.py -h` for info

## Visualization
### paraview.py:
To look at the density profiles in Sandia’s code Paraview, you can use the attached  script, e.g.:

- `python paraview Nx Ny Nz rhoga.tec rhoga.all`
-- From Jason Koski

-- Notes
Upload tec file
press apply
color -> real
put slice to see a plane
iso -> contour
