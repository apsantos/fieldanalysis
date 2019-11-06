#!/usr/bin/python
import math as m
import os, sys
import numpy as np

if ( len( sys.argv ) < 2 ):
  print "Usage: extra density [input.dat] [output.dat]"

else:

  inp = open(sys.argv[1] , 'r' )

  ndat = 0
  line = inp.readline().split()
  while ( len( line ) > 0):
    ndat += 1
    line = inp.readline().split()
    
  inp.close()

  print ndat
  
  inp = open( sys.argv[1] , 'r' )
  line = inp.readline().split()

  density = np.zeros((ndat,2),'d')
  for i in range( 0 , ndat ):
    density[i,0] = float(line[0])
    density[i,1] = float(line[1])
    line = inp.readline().split()
  den_sum = 0.0

  for i in range( 0 , ndat ):
    den_sum += density[i,1]

  print den_sum
  error = 100000
  tol = 1e-6
  while (abs(error) > tol):
    COM = 0.0
    for i in range(0,ndat):
      COM = COM + density[i,1] * (density[i,0] - m.floor(density[i,0]))
    COM = COM / den_sum
    error = COM - 0.5
    
    for i in range(0,ndat):
      density[i,0] = density[i,0] - error - m.floor(density[i,0] - error)
#  for i in range(0,ndat):
#    density[i,0] += 0.5
#    if (density[i,0] > 1.0):
#      density[i,0] -= 1.0

  density = density[np.argsort(density[:,0])] 
  otp = open(sys.argv[2],'w')
  for i in range(0,ndat):
    line = '%lf %lf\n' % (density[i,0], density[i,1])
    otp.write(line)

  otp.close()

