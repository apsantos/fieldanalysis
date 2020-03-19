#!/usr/bin/python
import sys, argparse
import math as ma
import numpy as np
import scipy
from scipy.special import erf
from scipy.optimize import curve_fit
import scipy.odr.odrpack as odrpack

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class FitProfile(object):
    #def erf_profile(self, x, rho_g, rho_l, xi, frac):
    def __init__(self):
        self.sqrt2 = ma.sqrt(2.0)

    def erf_profile_odr(self, B, x):
        rho_g = B[0]
        rho_l = B[1]
        xi = B[2]
        frac = B[3]
        factor = 0.5*(rho_l - rho_g) 
        return rho_g + factor * (erf((x - self.L_2 + frac) / xi / self.sqrt2) - erf((x - self.L_2 - frac) / xi / self.sqrt2))
        #return rho_g + 0.5*(rho_l - rho_g) * erf(ma.sqrt(ma.pi) * (x - self.L/2.0 + frac) / xi) - 0.5*(rho_l - rho_g) * erf(ma.sqrt(ma.pi) * (x - self.L/2.0 - frac) / xi)

    def erf_profile_curve(self, x, rho_g, rho_l, xi, frac):
        factor = 0.5*(rho_l - rho_g) 
        return rho_g + factor * (erf((x - self.L_2 + frac) / xi / self.sqrt2) - erf((x - self.L_2 - frac) / xi / self.sqrt2))

    def erf_set(self, x):
        factor = 0.5*(self.rho_l - self.rho_g) 
        return self.rho_g + factor * (erf((x - self.L_2 + frac) / self.xi / self.sqrt2) - erf((x - self.L_2 - frac) / self.xi / self.sqrt2))

    def tanh_profile(self, x, rho_g, rho_l, xi, frac):
        return 0.5*( rho_l + rho_g - (rho_l - rho_g) * np.tanh((x - self.L_2 - frac)/xi))
    
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Average meshed-field density'
                                                 'and convert to 1-d density profile')
    parser.add_argument("--profile_file", type=str, required=True, help='input density profile file')
    parser.add_argument("--fit_file", type=str, default='fitted_profile.dat', help='output fitted density profile file')
    parser.add_argument("--L_direction", type=float, default=1, help='length of simulation box in averaging direction')
    parser.add_argument('-l', "--start_line", type=int, default=1, help='starting line of density profile file')
    parser.add_argument("--fit_method", type=str, default='odr', choices=['odr','ls'], 
        help='odr: orthogonal distance regression, ls: least squares')

    L = parser.parse_args().L_direction

    ndat = file_len(parser.parse_args().profile_file) - parser.parse_args().start_line
    x = np.zeros(ndat,'d')
    y = np.zeros(ndat,'d')
#
    # read in position and density
    inp = open(parser.parse_args().profile_file, 'r') 

    i = 0
    j = 0
    i_half_a = 0
    i_half_b = 0
    for line in inp:
        data = line.strip().split()
        if j < parser.parse_args().start_line:
            j += 1
            continue
        x[i] = float(data[0])
        y[i] = float(data[1]) 
        if y[i] > 0.5 and i_half_a == 0:
            i_half_a = i
        if y[i] < 0.5 and i_half_b == 0:
            if i_half_a != 0:
                i_half_b = i
        i += 1

    guess= np.zeros(4,'d')
    guess[0] = y[0]
    guess[1] = y[int(ndat/2)]
    guess[2] = 0.1
    guess[3] = (x[i_half_b]-x[i_half_a])/2

    profile = FitProfile()
    profile.L = L # set the L in the function
    profile.L_2 = L/2.0 # set the L in the function

    fit_odr = False
    fit_ls = False
    if parser.parse_args().fit_method == 'odr': 
        fit_odr = True
    elif parser.parse_args().fit_method == 'ls': 
        fit_ls = True

    if fit_odr:
        model = odrpack.Model(profile.erf_profile_odr)
        data = odrpack.RealData(x,y)
        myodr = odrpack.ODR(data, model, beta0=guess,maxit=100,sstol=1e-5)
        output = myodr.run()
        fit_data = output.beta
        fit_stdev = output.sd_beta
 
    # fit the curve
    elif fit_ls:
        n = 0
        popt, pcov = curve_fit(profile.erf_profile_curve,x,y,guess,maxfev = 100000)
        while(n < 100):
            popt, pcov = curve_fit(profile.erf_profile_curve,x,y,popt,maxfev = 100000)
            n += 1
        fit_stdev = np.sqrt(np.diag(pcov))
        fit_data = popt

    profile.rho_g = fit_data[0]
    profile.rho_l = fit_data[1]
    profile.xi = fit_data[2]
    profile.frac = fit_data[3]
    print profile.rho_l
    print "# rho_g rho_l xi frac stdev[rho_g rho_l xi frac]"
    print fit_data[0], fit_data[1], fit_data[2], fit_data[3], fit_stdev[0], fit_stdev[1], fit_stdev[2], fit_stdev[3]

    fit = np.zeros(ndat,'d')
    for i in range(0,ndat):
        if fit_odr:
            fit[i] = profile.erf_profile_odr(fit_data, x[i])
        elif fit_ls:
            fit[i] = profile.erf_profile_curve(x[i], fit_data[0], fit_data[1], fit_data[2], fit_data[3])

    otp = open(parser.parse_args().fit_file, 'w') 
    otp.write('# position density\n')
    for i in range(0,ndat):
        line = '%lf %lf\n' % (x[i], fit[i])
        otp.write(line)

    otp.close()

if __name__ == '__main__':
    sys.exit(main())
