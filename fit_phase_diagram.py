#!/usr/local/bin/python
# !/usr/bin/python
import sys, argparse
import math as ma
import numpy as np
#import scipy
#from scipy.optimize import curve_fit
#import scipy.odr.odrpack as odrpack
from lmfit import minimize, Parameters, report_fit, Minimizer

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class FitPhaseCoex(object):
    def __init__(self):
        self.beta = 0.325
        self.mu =0.872
        #self.mu = 1.0 #0.872
        self.sign = 1.0

    def objective(self, params, x, data):
        """ calculate total residual for fits to several data sets held
        in a 2-D array, and modeled by Gaussian functions"""
        # make residual per data set
        ndata, nx = data.shape
        resid = 0.0*data[:]
        resid[0, :] = data[0, :] - self.thermo(params, x)
        resid[1, :] = data[1, :] - self.density(params, x)
        # now flatten this to a 1D array, as minimize() needs
        return resid.flatten()

    def density(self, B, T):
        rho_c = B['rho_c'].value
        T_c = B['T_c'].value
        delta_rho_0 = B['delta_rho_0'].value

        return delta_rho_0 * (self.sign*(T_c-T))**self.beta
        #return delta_rho_0 * ((T_c - T))**self.beta

    def thermo(self, B, T):
        rho_c = B['rho_c'].value
        T_c = B['T_c'].value
        A = B['A'].value

        return rho_c + (A * (self.sign*(T_c - T))**self.mu)

    def density_set(self, T):
        avg_rho = self.rho_c + self.A * self.sign * (self.T_c - T)
        delta_rho = self.delta_rho_0 * (self.sign * (self.T_c - T))**self.beta
        #delta_rho = self.delta_rho_0 * (self.sign * (1-(T/self.T_c)))**self.beta
        rho_l = avg_rho + delta_rho/2.0
        rho_v = avg_rho - delta_rho/2.0
        return rho_v, rho_l

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit vapor, liquid density to get critical point and curve')
    parser.add_argument("--profile_file", type=str, required=True, help='input density profile file')
    parser.add_argument("--fit_file", type=str, default='fitted_profile.dat', help='output fitted density profile file')
    parser.add_argument('-l', "--start_line", type=int, default=1, help='starting line of density profile file')
    parser.add_argument("--temp_col", type=int, default=0, help='temperature col')
    parser.add_argument("--vapor_col", type=int, default=1, help='vapor col')
    parser.add_argument("--liquid_col", type=int, default=2, help='liquid col')
    parser.add_argument("--type", type=str, choices=['lcst','ucst'],default='ucst', help='Is this an upper or lower critcal solution temperature transition?')
    parser.add_argument("--Tmin", type=float, help='Default is to ignore and use full T range.')
    parser.add_argument("--Tmax", type=float, help='Default is to ignore and use full T range.')

    c_t = parser.parse_args().temp_col
    c_v = parser.parse_args().vapor_col
    c_l = parser.parse_args().liquid_col
    ndat = file_len(parser.parse_args().profile_file) - parser.parse_args().start_line

    if parser.parse_args().Tmin:
        Tmin = parser.parse_args().Tmin
    if parser.parse_args().Tmax:
        Tmax = parser.parse_args().Tmax

    T = np.zeros(ndat,'d')
    rho_v = np.zeros(ndat,'d')
    rho_l = np.zeros(ndat,'d')
    fitT = []
    fitrho_v = []
    fitrho_l = []

    # read in temeprature and vapor/liquid densities
    inp = open(parser.parse_args().profile_file, 'r') 

    j = 0
    i = 0
    for line in inp:
        data = line.strip().split()
        if j < parser.parse_args().start_line:
            j += 1
            continue
        T[i] = float(data[c_t])
        rho_v[i] = float(data[c_v]) 
        rho_l[i] = float(data[c_l]) 
        if Tmin < T[i] and T[i] < Tmax:
            fitT.append( float(data[c_t]) )
            fitrho_v.append( float(data[c_v]) ) 
            fitrho_l.append( float(data[c_l]) )
        i += 1
    nfitdata = len(fitT)
    fitT = np.array(fitT)
    print fitT
    fit_data = np.zeros((2,nfitdata),'d')
    for i in range(nfitdata):
        rho_mean = (rho_v[i]+rho_l[i])/2.0
        delta_rho = rho_l[i] - rho_v[i]
        fit_data[:,i] = [rho_mean,delta_rho]

    if parser.parse_args().type == 'ucst':
        T_c_guess = max(T)*1.01
        T_c_min = T_c_guess
        if T_c_guess > 0:
            T_c_max = 5.0*T_c_guess
        else:
            T_c_max = 0.5*T_c_guess
    elif parser.parse_args().type == 'lcst':
        T_c_guess = min(T) * 0.99
        T_c_max = T_c_guess
        if T_c_guess > 0:
            T_c_min = 0.5*T_c_guess
        else:
            T_c_min = 5.0*T_c_guess

    rho_c_guess = (max(rho_v)+min(rho_l))/2.0

    print T_c_min, T_c_max
    fit_params = Parameters()
    fit_params.add( 'rho_c', value=rho_c_guess, min=max(rho_v),  max=min(rho_l))
    fit_params.add( 'T_c', value=T_c_guess, min=T_c_min, max=T_c_max)
    fit_params.add( 'A', value=0.1, min=0.0,  max=100.0)
    fit_params.add( 'delta_rho_0', value=0.1, min=0.0, max=3.0)
    
    profile = FitPhaseCoex()
    if parser.parse_args().type == 'ucst':
        profile.sign = 1.0
    elif parser.parse_args().type == 'lcst':
        profile.sign = -1.0
    
    minner = Minimizer(profile.objective, fit_params, fcn_args=(fitT, fit_data), calc_covar=True)
    #minner = Minimizer(profile.objective, fit_params, fcn_args=(T, fit_data), calc_covar=True)
    result = minner.minimize()
    #report_fit(result)

    profile.rho_c = result.params['rho_c'].value
    profile.T_c = result.params['T_c'].value
    profile.A = result.params['A'].value
    profile.delta_rho_0 = result.params['delta_rho_0'].value

    # calculate final result
    # write critical point
    otp = open(parser.parse_args().fit_file+".crit", 'w') 
    otp.write('# \\frac{\\rho_v+\\rho_l}{2} = \\rho_c + A(1-T/T_c)\n')
    otp.write('# \\delta\\rho = \\delta\\rho_0(1-T/T_c)^{\\beta}\n')
    otp.write('# Fit quality chi^2: %lf chi^2/N: %lf\n' % (result.chisqr, result.redchi) )
    otp.write('# T_c rho_c A \delta\\rho_0')
    if hasattr(result, 'covar'): otp.write(' stdev[T_c rho_c A \delta\\rho_0]')
    otp.write('\n')
    otp.write('%lf %lf %lf %lf' % (profile.T_c, profile.rho_c, profile.A, profile.delta_rho_0) )
    if hasattr(result, 'covar'):
        otp.write(' %lf %lf %lf %lf' % (result.covar[1,1], result.covar[0,0], result.covar[2,2], result.covar[3,3]) )
    otp.write('\n')
    otp.close()

    # write density fit
    nfit = 500
    if parser.parse_args().type == 'ucst': 
        Tlow = min(T)
        Tlow = Tmin
    elif parser.parse_args().type == 'lcst': 
        Tlow = max(T)
        Tlow = Tmax
    fit_T = np.linspace(Tlow,profile.T_c,num=nfit)

    otp = open(parser.parse_args().fit_file, 'w') 
    otp.write('# T rho_v rho_l\n') 
    for i in range(nfit):
        fit_v, fit_l = profile.density_set(fit_T[i])
        line = '%lf %lf %lf\n' % (fit_T[i], fit_v, fit_l)
        otp.write(line)

    otp.close()

if __name__ == '__main__':
    sys.exit(main())

