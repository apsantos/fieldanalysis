#!/usr/bin/env python3
# #!/Users/asanto/anaconda3/bin/python
# #!/usr/bin/python
import sys, argparse
import math as ma
import numpy as np
import scipy
from scipy.special import erf
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odrpack
from ovito.io import import_file, export_file
from ovito.modifiers import CombineDatasetsModifier, DeleteSelectedModifier, AffineTransformationModifier, ExpressionSelectionModifier, InvertSelectionModifier

# User-defined modifier 'Shrink-wrap simulation box':
def modify(frame, data):

    # There's nothing we can do if there are no particles. 
    if data.particles.count == 0: return

    # Compute bounding box of particle coordinates.
    coords_min = np.amin(data.particles.positions, axis=0)
    coords_max = np.amax(data.particles.positions, axis=0)

    # Adjust simulation cell matrix (three cell vectors and origin).
    data.cell_[:,:3] = np.diag(coords_max - coords_min)
    data.cell_[:, 3] = coords_min
# This is a copy of the template file '/Applications/Ovito.app/Contents/Resources/scripts/modifiers/Shrink-wrap simulation box.py'.
# Feel free to modify the code below as needed.

#
# Shrink-wrap simulation box:
#
# A user-defined modifier function which resets the simulation box geometry to 
# exactly match the current axis-aligned bounding box of the particle coordinates.
#
   
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Average meshed-field density'
                                                 'and convert to 1-d density profile')
    parser.add_argument("--configfile", type=str, required=True, help='configuration filename')
    parser.add_argument("--outputfile", type=str, required=True, help='filename with crystallinity as a function of time')
    parser.add_argument("--left", type=int, required=True, help='left cut')
    parser.add_argument("--right", type=int, required=True, help='right cut')

    # Data import:
    ifile =  parser.parse_args().configfile
    ofile =  parser.parse_args().outputfile
    left =  parser.parse_args().left
    right =  parser.parse_args().right
    shift = right - left

    pipeliner = import_file(ifile, atom_style = 'molecular')
    pipeliner.modifiers.append(ExpressionSelectionModifier(expression = 'Position.Z>'+str(left)))
    pipeliner.modifiers.append(InvertSelectionModifier())
    pipeliner.modifiers.append(DeleteSelectedModifier())
    pipeliner.modifiers.append(AffineTransformationModifier(transformation = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, shift]]))
    export_file(pipeliner, "right.data", "lammps/data", atom_style = "molecular")

    pipelinel = import_file(ifile, atom_style = 'molecular')
    pipelinel.modifiers.append(ExpressionSelectionModifier(expression = 'Position.Z<'+str(right)))
    pipelinel.modifiers.append(InvertSelectionModifier())
    pipelinel.modifiers.append(DeleteSelectedModifier())
    export_file(pipelinel, "left.data", "lammps/data", atom_style = "molecular")

    # read in left side
    pipelinen =  import_file('left.data', atom_style = 'molecular')
    # combine right side
    mod = CombineDatasetsModifier()
    mod.source.load('right.data', atom_style = 'molecular')
    pipelinen.modifiers.append(mod)
    # shrink wrap box
    pipelinen.modifiers.append(modify)
    export_file(pipelinen, ofile, "lammps/data", atom_style = "molecular")

if __name__ == '__main__':
    sys.exit(main())
