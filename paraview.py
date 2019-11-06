#!/usr/bin/python

import math, sys, string, copy, random, os
from sys import stderr, stdout, stdin
import numpy as nmp
from math import *

def main():

  if ( len( sys.argv ) < 5 ):
    print ">>>> These are the inputs [Nx], [Ny], [Nz], [rho.tec] [rho.all]"
    exit(1)

  Nx = int(sys.argv[1])
  Ny = int(sys.argv[2])
  Nz = int(sys.argv[3])
  File = sys.argv[4]

  hdr = open(sys.argv[4], 'w')
  hdr.write("TITLE = \"%s\"\n" % File)
  hdr.write("VARIABLES = \"X\", \"Y\", \"Z\", \"Real\"" ); 
  hdr.write("ZONE I=%d, J=%d, K=%d, F=POINT\n" % (Nx, Ny, Nz) ); 
  hdr.close()

  cmd = 'cat %s >> %s' % ( sys.argv[5] , sys.argv[4] )
  os.system(cmd)

  return

if __name__ == "__main__":
  main()

