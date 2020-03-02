#!/usr/bin/python
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

def func_fit(x,rho_g, rho_l, xi, frac):
  #print rho_g, rho_l, xi, frac  
  return rho_g + 0.5*(rho_l - rho_g) * scipy.special.erf(m.sqrt(m.pi) * (x - Lx/2.0 + frac) / xi) - 0.5*(rho_l - rho_g) * scipy.special.erf(m.sqrt(m.pi) * (x - Lx/2.0 - frac) / xi)

if ( len( sys.argv ) < 3 ):
    print "Usage: radial_avg_rho.py [Lx] [input.dat] [output.dat] [iter number]"

else:
  Lx = float(sys.argv[1])
  inp = open(sys.argv[2], 'r') 

  line = inp.readline().split()

  ndat = 0

  while (len(line) > 0):
    if float(line[1]) > 0.0:
      ndat += 1
    line = inp.readline().split()

  #print ndat

  inp.close()

  x = np.zeros(ndat,'d')
  y = np.zeros(ndat,'d')

  inp = open(sys.argv[2], 'r')

  line = inp.readline().split()

  for i in range(0,ndat):
    x[i] = float(line[0])
    y[i] = float(line[1]) 
    line = inp.readline().split()

  guess= np.zeros(4,'d')
  guess[0] = 0.1
  guess[1] = 0.9
  guess[2] = 0.1
  guess[3] = 0.1
  #print guess

  popt, pcov = curve_fit(func_fit,x,y,guess,maxfev = 100000)
 
  #print popt;

  n = 0
  while(n < 100):
    popt, pcov = curve_fit(func_fit,x,y,popt,maxfev = 100000)
    n += 1
  print popt[0], popt[1], popt[2], popt[3], pcov[0], pcov[1], pcov[2], pcov[3]

  fit = np.zeros(ndat,'d')
  for i in range(0,ndat):
    fit[i] = popt[0] + 0.5*(popt[1] - popt[0]) * scipy.special.erf(m.sqrt(m.pi) * (x[i] - Lx/2.0 + popt[3]) / popt[2]) - 0.5*(popt[1] - popt[0]) * scipy.special.erf(m.sqrt(m.pi) * (x[i] - Lx/2.0 - popt[3]) / popt[2])

  otp = open(sys.argv[3], 'w')
  for i in range(0,ndat):
    line = '%lf %lf\n' % (x[i], fit[i])
    otp.write(line)

  iter_num = int(sys.argv[4])
  #print 100*100*(popt[2] * popt[2] / 2.0 / m.pi)
  #print popt[0]
  #print popt[1]


  otp.close()
