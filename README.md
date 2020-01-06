# From Jason Koski

## Calculate the density profiles in the z direction and fit an error function to the two interfaces:
### dmft_slice.py: 
Python script that takes an output density file and averages over a chosen plane.
- dmft_slice.py Nx Ny Nz L_{chosen_plane} {chosen_plane [0:x, 1:y, 2:z]} {file doing analysis on} {rho_0}
-- Nx = Number of grid points in x-direction
-- Lx = Length of box in x-direction
-- Hardcoded output is always 'dens_no-shift_slice.dat'

### com.py: 
Shifts input file so that the density file is centered at Lz/2

- com.py [input file, typically 'dens_no-shift_slice.dat'] [centered output density]

- com2.py is basically the same file as com.py but it shifts the density file to be centered at 0

Once in a while in doing my mass analysis, some of the density profiles would be shifted so that the matrix-dense region was more dense in the center of the box. 'fit_curve.py' is best setup to calculate density profiles for [matrix-dense, gnp-dense, matrix-dense] curves. This way, you have both centering scripts incase you run into this problem.

- fit_curve.py: fits a density curve to 2 error functions and spits out the matrix-dense volume fraction and the gnp-dense volume fraction
-- fit_curve.py Lz [input density curve] [output fitted curve] [dummy variable, just set this to 1, or whatever, its not used]
--Note dmft_slice.py spits out the density file as phi vs z/Lz, so always use Lz is 1
--If you get funky results from this, you likely have to center the density curve using com2.py OR once in a blue moon, you may have to change your initial guess.
--Note you’ll need the scripts in the fts folder as well, as they have some functions that are used.

### paraview.py:
To look at the density profiles in Sandia’s code Paraview, you can use the attached  script, e.g.:

- python paraview Nx Ny Nz rhoga.tec rhoga.all

Notes
Upload tec file
press apply
color -> real
put slice to see a plane
iso -> contour
