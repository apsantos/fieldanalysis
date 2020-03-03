#!/usr/bin/python
import sys, argparse
import math as m
import cmath as cm
import os, sys
import numpy as np
import scipy
import random
from scipy.special import erf
from scipy.fftpack import fft, ifft
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class ERFprofile(object):
    def erf_profile(self, x, rho_g, rho_l, xi, frac):
        return rho_g + 0.5*(rho_l - rho_g) * scipy.special.erf(m.sqrt(m.pi) * (x - self.L/2.0 + frac) / xi) - 0.5*(rho_l - rho_g) * scipy.special.erf(m.sqrt(m.pi) * (x - self.L/2.0 - frac) / xi)
    
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Average meshed-field density'
                                                 'and convert to 1-d density profile')
    parser.add_argument("--profile_file", type=str, required=True, help='input density profile file')
    parser.add_argument("--fit_file", type=str, default='fitted_profile.dat', help='output fitted density profile file')
    parser.add_argument("--L_direction", type=float, default=1, help='length of simulation box in averaging direction')

    L = parser.parse_args().L_direction

    ndat = file_len(parser.parse_args().profile_file)

    x = np.zeros(ndat,'d')
    y = np.zeros(ndat,'d')

    inp = open(parser.parse_args().profile_file, 'r') 

    for i in range(0,ndat):
      line = inp.readline().split()
      x[i] = float(line[0])
      y[i] = float(line[1]) 

    guess= np.zeros(4,'d')
    guess[0] = 0.1
    guess[1] = 0.9
    guess[2] = 0.1
    guess[3] = 0.1

    # set the L in the function
    profile = ERFprofile()
    profile.L = L
    popt, pcov = curve_fit(profile.erf_profile,x,y,guess,maxfev = 100000)
 
    n = 0
    while(n < 100):
      popt, pcov = curve_fit(profile.erf_profile,x,y,popt,maxfev = 100000)
      n += 1
    print "# rho_g rho_l xi frac uncertainties
    print popt[0], popt[1], popt[2], popt[3], pcov[0], pcov[1], pcov[2], pcov[3]

    fit = np.zeros(ndat,'d')
    for i in range(0,ndat):
      #fit[i] = popt[0] + 0.5*(popt[1] - popt[0]) * scipy.special.erf(m.sqrt(m.pi) * (x[i] - L/2.0 + popt[3]) / popt[2]) - 0.5*(popt[1] - popt[0]) * scipy.special.erf(m.sqrt(m.pi) * (x[i] - L/2.0 - popt[3]) / popt[2])
      fit[i] = profile.erf_profile(x[i], popt[0], popt[1], popt[2], popt[3])

    otp = open(parser.parse_args().fit_file, 'w') 
    for i in range(0,ndat):
      line = '%lf %lf\n' % (x[i], fit[i])
      otp.write(line)

    otp.close()

if __name__ == '__main__':
    sys.exit(main())
