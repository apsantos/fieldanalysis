#!/usr/bin/python
import sys, argparse
import math as ma
import numpy as np


def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='NP-NP separation prediction')
    parser.add_argument("--outfile", type=str, default='prediction.rnpnp', help='output fitted density profile file')
    parser.add_argument("--Dnp", type=float, required=True, help='NP diameter')
    parser.add_argument("--maxphi", type=float, required=True, default=0.64, help='maximum phi (eg RCP:0.64, HCP:0.74, BCC:0.68')
    parser.add_argument("--norm", type=float, default=1.0, help='normalize distance')
    parser.add_argument("--pdi", type=float, default=1.0, help='breadth  of  the  log-normal size distribution (eg 1 for monodisperse')

    out = open(parser.parse_args().outfile, 'w')
    Dnp = parser.parse_args().Dnp
    maxphi = parser.parse_args().maxphi
    pdi = parser.parse_args().pdi
    norm = parser.parse_args().norm

    j = 0
    factor1 = Dnp*ma.exp(1.5*np.log(pdi)**2.)
    factor2 = 0.0 #Dnp*ma.exp(0.5*np.log(pdi)**2.)
    print (maxphi/0.0954)**(1./3.) * factor1 - factor2, factor1, factor2, (maxphi/0.0954)**(1./3.), maxphi
    for phi in np.arange(0.001,0.25,0.001):
        rnpnp = (maxphi/phi)**(1./3.) * factor1 - factor2
        rnpnp_norm = rnpnp / norm
        out.write("%e %e %e\n" % (phi, rnpnp, rnpnp_norm) )

    out.close()
if __name__ == '__main__':
    sys.exit(main())
