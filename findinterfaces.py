#!/usr/bin/python
# #!/Users/asanto/anaconda3/bin/python
# #!/usr/bin/env python3
import sys, argparse
import math as ma
import numpy as np
from scipy.special import erf
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar

def erffit(z, rho_goal, rho_g, rho_l, xi, frac):
    xisqrt = xi * ma.sqrt(2)
    rhoi = rho_g + (0.5 * (rho_l - rho_g) * (erf((z - 0.5 + frac)/xisqrt) - erf((z - 0.5 - frac)/xisqrt)))
    return rhoi - rho_goal

class erffitz(object):
    def func(self,z):
        rhoi = self.rho_g + (0.5 * (self.rho_l - self.rho_g) * (erf((z - 0.5 + self.frac)/self.xisqrt) - erf((z - 0.5 - self.frac)/self.xisqrt))) - self.rho_goal
        return abs(rhoi)
    
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Average meshed-field density'
                                                 'and convert to 1-d density profile')
    parser.add_argument("--interfacecutoff", type=float, default=0.02, 
                            help='percentage of yfit_low and yfit_hi to start interface')
    parser.add_argument("--yfit_low", type=float, help='erf fit low phi/y value')
    parser.add_argument("--yfit_hi", type=float, help='erf fit low phi/y value')
    parser.add_argument("--xifit", type=float, help='minimum bound for the hi phi fit value')
    parser.add_argument("--fracfit", type=float, help='minimum bound for the low phi fit value')

    interfacecutoff = parser.parse_args().interfacecutoff
    yfit_low = parser.parse_args().yfit_low
    yfit_hi = parser.parse_args().yfit_hi
    xifit = parser.parse_args().xifit
    fracfit = parser.parse_args().fracfit

    ef = erffitz()
    ef.rho_g = yfit_low
    ef.rho_l = yfit_hi
    ef.frac = fracfit
    ef.xisqrt = xifit * ma.sqrt(2)

    delta = yfit_hi - yfit_low
    dcut = delta * interfacecutoff

    ef.rho_goal = yfit_low + dcut
    result = minimize_scalar( ef.func, bounds=(0.0,0.5), method='bounded' )
    zla = result.x

    ef.rho_goal = yfit_hi - dcut
    result = minimize_scalar( ef.func, bounds=(zla,0.5), method='bounded' )
    zha = result.x
    print zla, zha, 1-zha, 1-zla

if __name__ == '__main__':
    sys.exit(main())
