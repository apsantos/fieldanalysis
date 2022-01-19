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
from ovito.modifiers import SelectTypeModifier, DeleteSelectedModifier, AcklandJonesModifier


   
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Average meshed-field density'
                                                 'and convert to 1-d density profile')
    parser.add_argument("--configfiles", type=str, required=True, nargs='+', help='configuration filename')
    parser.add_argument("--outputfile", type=str, required=True, help='filename with crystallinity as a function of time')
    parser.add_argument("--timestep", type=float, default=0.02, help='timestep')

    pipeline = import_file(parser.parse_args().configfiles)
    # Select type:
    pipeline.modifiers.append(SelectTypeModifier(types = {10, 4, 5, 15}))
    
    # Delete selected:
    pipeline.modifiers.append(DeleteSelectedModifier())
    
    # Ackland-Jones analysis:
    pipeline.modifiers.append(AcklandJonesModifier())
    export_file(pipeline, parser.parse_args().outputfile, "txt/attr", multiple_frames = True,
    columns = ["Timestep", 'AcklandJones.counts.FCC','AcklandJones.counts.HCP','AcklandJones.counts.BCC','AcklandJones.counts.ICO','AcklandJones.counts.OTHER'] )

    # calculte precentages
    ifile = open(parser.parse_args().outputfile, 'r' )
    ofile = open(parser.parse_args().outputfile + 'norm', 'w' )
    line = ifile.readline()
    ofile.write( '# time % FCC HCP BCC ICO OTHER\n' )
    pstep = 0
    prstep = 0
    for line in ifile:
        data = line.strip().split()
        step = int(data[0]) 
        if step == 0: prstep += pstep
        time = (step + prstep) * parser.parse_args().timestep
        fcc = int(data[1])
        hcp = int(data[2])
        bcc = int(data[3])
        ico = int(data[4])
        other = int(data[5])
        total = float( fcc + hcp + bcc + ico + other )
        fcc /= total
        hcp /= total
        bcc /= total
        ico /= total
        other /= total
        ofile.write( '%f %f %f %f %f %f\n' % (time, fcc, hcp, bcc, ico, other) )
        pstep = step

    ofile.close()
    ifile.close()

if __name__ == '__main__':
    sys.exit(main())
