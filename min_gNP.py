#!/usr/bin/python
import sys, argparse
from scipy.special import erfinv
import math as ma

def calcPhaseVolume(box, direction, width, z0, fraction=0.9975):
    if   direction == 'x': dim = 0
    elif direction == 'y': dim = 1
    elif direction == 'z': dim = 2

    # calculate gNP-poor phase volume
    vol = 1.0
    Ldilute = 2*( (-width*ma.sqrt(2)*erfinv(fraction)) + 0.5 - z0)
    for i in [0,1,2]:
        if i == dim:
            vol *= Ldilute*box[i]
        else:
            vol *= box[i]
    return vol
        
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Read in mesh/solvent 1-d density profile and'
                                                 'determine if there are enough gNP for accurate gNP-poor region measurement')
    parser.add_argument("--direction", type=str, choices=['x', 'y', 'z'], default='z', help='1-d profile direction')
    parser.add_argument("--L", type=float, nargs='+', help='length of sim box, eg: --L 25 25 50')
    parser.add_argument("--fitted_phiM", type=float, required=True, help='gNP-poor volume fraction')
    parser.add_argument("--fitted_width", type=float, required=True, help='fit width')
    parser.add_argument("--fitted_z0", type=float, required=True, help='fit width')
    parser.add_argument("--NgNP", type=int, required=True, help='number of grafted nanoparticles')
    parser.add_argument("--Rp", type=float, required=True, help='NP radius [b]')
    parser.add_argument("--sigma", type=float, required=True, help='graft density [chains/b^2]')
    parser.add_argument("--NG", type=int, required=True, default=1, help='graft length')

    NgNP = parser.parse_args().NgNP
    Rp = parser.parse_args().Rp
    sigma = parser.parse_args().sigma
    NG = parser.parse_args().NG
    fit_phiM = parser.parse_args().fitted_phiM

    vol = calcPhaseVolume(parser.parse_args().L, parser.parse_args().direction, parser.parse_args().fitted_width, parser.parse_args().fitted_z0)
    min_phi_gNP = ((4.*ma.pi / 3. * Rp**3.) + (ma.sqrt(NG/6) * sigma * 4.*ma.pi * Rp**2.)) / vol
    print 'min(phi_gNP), phi_gNP:', min_phi_gNP, (1.-fit_phiM)

if __name__ == '__main__':
    sys.exit(main())
