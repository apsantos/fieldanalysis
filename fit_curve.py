#!/Users/asanto/anaconda3/bin/python
# #!/usr/bin/python
# #!/usr/bin/env python3
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
        rho_g = max(B[0],self.minB)
        rho_l = min(B[1],self.maxB)
        for i in [2,3]:
            B[i] = max(B[i],0.00001)
            B[i] = min(B[i],1)
        xi = B[2]
        frac = B[3]
        factor = 0.5*(rho_l - rho_g) 
        factor2 = xi * self.sqrt2
        return rho_g + factor * (erf((x - self.L_2 + frac) / factor2) - erf((x - self.L_2 - frac) / factor2))
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
    parser.add_argument("--fit_method", type=str, default='odr', choices=['odr','trf','lm'], 
        help='odr: orthogonal distance regression, lm: least squares')
    parser.add_argument("--phi_low_max", type=float, help='maximum bound for the low phi fit value')
    parser.add_argument("--phi_hi_min", type=float, help='minimum bound for the hi phi fit value')
    parser.add_argument("--phi_low_min", type=float, help='minimum bound for the low phi fit value')
    parser.add_argument("--phi_hi_max", type=float, help='maximum bound for the hi phi fit value')

    L = parser.parse_args().L_direction

    set_phi_low_min = False
    if parser.parse_args().phi_low_min:
        set_phi_low_min = True
        phi_low_min = parser.parse_args().phi_low_min

    set_phi_low_max = False
    if parser.parse_args().phi_low_max:
        set_phi_low_max = True
        phi_low_max = parser.parse_args().phi_low_max

    set_phi_hi_min = False
    if parser.parse_args().phi_hi_min:
        set_phi_hi_min = True
        phi_hi_min = parser.parse_args().phi_hi_min

    set_phi_hi_max = False
    if parser.parse_args().phi_hi_max:
        set_phi_hi_max = True
        phi_hi_max = parser.parse_args().phi_hi_max

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
        i += 1

    ymin = min(y)
    ymax = max(y)
    #ymin = 1.-max(y)
    #ymax = 1.-min(y)
    halfway = ymin + 0.5*(ymax - ymin)
    for i in range(len(y)):
        if y[i] > halfway and i_half_a == 0:
            i_half_a = i
            break
    for i in range(i_half_a+20,len(y)):
        if y[i] < halfway and i_half_b == 0:
            i_half_b = i
            break

    guess= np.zeros(4,'d')
    guess[0] = max(ymin,0.0)
    guess[1] = min(ymax,1.0)
    guess[2] = 0.5
    guess[3] = (x[i_half_b]-x[i_half_a])/2.

    profile = FitProfile()
    profile.L = L # set the L in the function
    profile.L_2 = L/2.0 # set the L in the function

    fit_odr = False
    fit_lm = False
    fit_trf = False
    if parser.parse_args().fit_method == 'odr': 
        fit_odr = True
    elif parser.parse_args().fit_method == 'lm': 
        fit_lm = True
    elif parser.parse_args().fit_method == 'trf': 
        fit_trf = True

    if fit_odr:
        profile.minB = 0.000001
        profile.maxB = 1.0
        if set_phi_low_min: profile.minB = phi_low_min
        if set_phi_hi_max:  profile.maxB = phi_hi_max
        
        model = odrpack.Model(profile.erf_profile_odr)
        data = odrpack.RealData(x,y)
        myodr = odrpack.ODR(data, model, beta0=guess,maxit=1000,sstol=1e-6)
        output = myodr.run()
        fit_data = output.beta
        fit_stdev = output.sd_beta
 
    # fit the curve
    elif fit_lm or fit_trf:
        if set_phi_low_min: 
            plowmin = phi_low_min
        else:
            plowmin = 0.0000
        if set_phi_low_max: 
            plowmax = phi_low_max
        else:
            plowmax = 1.0000
        if set_phi_hi_min:  
            phimin = phi_hi_min
        else:
            phimin = 0.0000
        if set_phi_hi_max:  
            phimax = phi_hi_max
        else:
            phimax = 1.0
        bounds=([plowmin,phimin,0.00001,0.00001], [plowmax,phimax, 1.,1.])
        print (bounds)
        #print (guess)
        n = 0
        popt, pcov = curve_fit(profile.erf_profile_curve, x, y, guess, bounds=bounds, maxfev = 100000, method=parser.parse_args().fit_method)
        while(n < 100):
            popt, pcov = curve_fit(profile.erf_profile_curve, x, y, guess, bounds=bounds, maxfev = 100000, method=parser.parse_args().fit_method)
            n += 1
        fit_stdev = np.sqrt(np.diag(pcov))
        fit_data = popt

    profile.rho_g = fit_data[0]
    profile.rho_l = fit_data[1]
    profile.xi = fit_data[2]
    profile.frac = fit_data[3]
    print ("# rho_g rho_l xi frac stdev[rho_g rho_l xi frac]")
    print (fit_data[0], fit_data[1], fit_data[2], fit_data[3], fit_stdev[0], fit_stdev[1], fit_stdev[2], fit_stdev[3])

    fit = np.zeros(ndat,'d')
    for i in range(0,ndat):
        if fit_odr:
            fit[i] = profile.erf_profile_odr(fit_data, x[i])
        elif fit_lm or fit_trf:
            fit[i] = profile.erf_profile_curve(x[i], fit_data[0], fit_data[1], fit_data[2], fit_data[3])

    otp = open(parser.parse_args().fit_file, 'w') 
    otp.write('# position density\n')
    for i in range(0,ndat):
        line = '%lf %lf\n' % (x[i], fit[i])
        otp.write(line)

    otp.close()

if __name__ == '__main__':
    sys.exit(main())
