import os,sys
import numpy as nmp
import math as m
import cmath as cm

def calc_debye(x):
  return ( (cm.exp(-x*x)+x*x-1.0)/(x*x*x*x) )

def mod_avg_3d_full_slice_no_shift(Nx, Ny, Nz, loop_dir, rho):

  Nmx = nmp.max([Nx, Ny , Nz])

  rh = nmp.zeros(Nmx,'cdouble')
  rhstd = nmp.zeros(Nmx,'cdouble')
  nrm = nmp.zeros(Nmx, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz


  for k in range(int(L1/3) , int(2*L1/3)):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      for j in range(0 , LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rh[ j ] += rho[ind]
 
        nrm[ j ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rh[i] = rh[i] / nrm[i]


  for k in range(int(L1/3) , int(2*L1/3)):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      for j in range(0 , LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rhstd[ j ] += (rh[ j ] - rho[ind])**2
 
        nrm[ j ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rhstd[i] = (rhstd[i] / nrm[i]) ** 0.5


  return rh , rhstd

def avg_3d_orien_full_slice_no_shift(Nx, Ny, Nz, loop_dir, sxx, syy, szz, sxy, sxz, syz):

  Nmx = nmp.max([Nx, Ny , Nz])

  rhxx = nmp.zeros(Nmx,'cdouble')
  rhyy = nmp.zeros(Nmx,'cdouble')
  rhzz = nmp.zeros(Nmx,'cdouble')
  rhxy = nmp.zeros(Nmx,'cdouble')
  rhxz = nmp.zeros(Nmx,'cdouble')
  rhyz = nmp.zeros(Nmx,'cdouble')
  #rhstd = nmp.zeros(Nmx,'cdouble')
  nrmxx = nmp.zeros(Nmx, 'd')
  nrmyy = nmp.zeros(Nmx, 'd')
  nrmzz = nmp.zeros(Nmx, 'd')
  nrmxy = nmp.zeros(Nmx, 'd')
  nrmxz = nmp.zeros(Nmx, 'd')
  nrmyz = nmp.zeros(Nmx, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      for j in range(0 , LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rhxx[ j ] += sxx[ind]
        rhyy[ j ] += syy[ind]
        rhzz[ j ] += szz[ind]
        rhxy[ j ] += sxy[ind]
        rhxz[ j ] += sxz[ind]
        rhyz[ j ] += syz[ind]
 
        nrmxx[ j ] += 1.0
        nrmyy[ j ] += 1.0
        nrmzz[ j ] += 1.0
        nrmxy[ j ] += 1.0
        nrmxz[ j ] += 1.0
        nrmyz[ j ] += 1.0


  for i in range(0, LL):
    if (nrmxx[i] > 0.0):
      rhxx[i] = rhxx[i] / nrmxx[i]
    if (nrmyy[i] > 0.0):
      rhyy[i] = rhyy[i] / nrmyy[i]
    if (nrmzz[i] > 0.0):
      rhzz[i] = rhzz[i] / nrmzz[i]
    if (nrmxy[i] > 0.0):
      rhxy[i] = rhxy[i] / nrmxy[i]
    if (nrmxz[i] > 0.0):
      rhyz[i] = rhxz[i] / nrmxz[i]
    if (nrmyz[i] > 0.0):
      rhyz[i] = rhyz[i] / nrmyz[i]


 # for k in range(0 , L1):

 #   for i in range(0, L2):
 #
 #     ind2 = k + i*L1
 #     
 #     for j in range(0 , LL):

 #       ind = 0
 #
 #       if (loop_dir==0):
 #         ind = j + (k + i*Ny)*Nx
 #       elif (loop_dir == 1):
 #         ind = k + (j + i*Ny)*Nx
 #       else:
 #         ind = k + (i + j*Ny)*Nx

 #       rhstd[ j ] += (rh[ j ] - rho[ind])**2
 #
 #       nrm[ j ] += 1.0


 # for i in range(0, LL):
 #   if (nrm[i] > 0.0):
 #     rhstd[i] = (rhstd[i] / nrm[i]) ** 0.5


  return rhxx, rhyy, rhzz, rhxy, rhxz, rhyz

def avg_3d_full_slice_no_shift(Nx, Ny, Nz, loop_dir, rho):

  Nmx = nmp.max([Nx, Ny , Nz])

  rh = nmp.zeros(Nmx,'cdouble')
  rhstd = nmp.zeros(Nmx,'cdouble')
  nrm = nmp.zeros(Nmx, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      for j in range(0 , LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rh[ j ] += rho[ind]
 
        nrm[ j ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rh[i] = rh[i] / nrm[i]


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      for j in range(0 , LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rhstd[ j ] += (rh[ j ] - rho[ind])**2
 
        nrm[ j ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rhstd[i] = (rhstd[i] / nrm[i]) ** 0.5


  return rh , rhstd


def avg_2d_full_slice_no_shift(Nx, Ny, loop_dir, rho):

  Nmx = nmp.max([Nx, Ny ])

  rh = nmp.zeros(Nmx,'cdouble')
  rhstd = nmp.zeros(Nmx,'cdouble')
  nrm = nmp.zeros(Nmx, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny


  for k in range(0 , L1):

    for j in range(0 , LL):

      ind = 0
 
      if (loop_dir==0):
        ind = j + k * Nx
      elif (loop_dir == 1):
        ind = k + j * Nx

      rh[ j ] += rho[ind]
 
      nrm[ j ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rh[i] = rh[i] / nrm[i]


  for k in range(0 , L1):

    for j in range(0 , LL):
      
      ind = 0
 
      if (loop_dir==0):
        ind = j + k * Nx
      elif (loop_dir == 1):
        ind = k + j * Nx

      
      rhstd[ j ] += (rh[ j ] - rho[ind])**2
 
      nrm[ j ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rhstd[i] = (rhstd[i] / nrm[i]) ** 0.5


  return rh , rhstd



def dump_3d_dens_slice(Nx, Ny, Nz, loop_dir, rho):

  rh = nmp.zeros(Nx,'cdouble')
  rhstd = nmp.zeros(Nx,'cdouble')
  nrm = nmp.zeros(Ny, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz

  otp = open("dens_slice_dump.dat", 'w')

  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1

      for j in range(0, LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        line = '%d %e %e\n' % (j , rho[ind].real , rho[ind].imag )
        otp.write(line)

      otp.write('\n')

  otp.close()



def avg_3d_full_slice(Nx, Ny, Nz, loop_dir, xbgn, xwid , rho):

  rh = nmp.zeros(Nx,'cdouble')
  rhstd = nmp.zeros(Nx,'cdouble')
  nrm = nmp.zeros(Ny, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz

  avg_shift = 0.0;
  for i in range( 0 , L1*L2 ):
    avg_shift += xbgn[i][1] - xbgn[i][0]

  avg_shift = int( avg_shift / L1 / L2 )


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      minim = xbgn[ind2][1] - avg_shift - 10
      if (minim < 0):
        minim = 0

      maxim = xbgn[ind2][1] + 10
      if (maxim >= LL):
        maxim = LL

      for j in range(minim, maxim):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rh[ j-minim ] += rho[ind]
 
        nrm[ j-minim ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rh[i] = rh[i] / nrm[i]


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      minim = xbgn[ind2][1] - avg_shift - 7
      if (minim < 0):
        minim = 0
      
      maxim = xbgn[ind2][1] + 7
      if (maxim >= LL):
        maxim = LL


      for j in range(minim, maxim):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rhstd[ j-minim ] += (rh[ j-minim ] - rho[ind])**2
 
        nrm[ j-minim ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rhstd[i] = (rhstd[i] / nrm[i]) ** 0.5


  return rh , rhstd


def avg_3d_dens_slice(Nx, Ny, Nz, loop_dir, xbgn, rho):

  rh = nmp.zeros(Nx,'cdouble')
  rhstd = nmp.zeros(Nx,'cdouble')
  nrm = nmp.zeros(Ny, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1

      for j in range(xbgn[ind2], LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rh[ j-xbgn[ind2] ] += rho[ind]
 
        nrm[ j-xbgn[ind2] ] += 1.0

  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rh[i] = rh[i] / nrm[i]



  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      for j in range(xbgn[ind2] , LL):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rhstd[ j-xbgn[ind2] ] += (rh[ j-xbgn[ind2] ] - rho[ind])**2
 
        nrm[ j-xbgn[ind2] ] += 1.0

  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rhstd[i] = (rhstd[i] / nrm[i]) ** 0.5

  return rh , rhstd


def avg_3d_full_slice(Nx, Ny, Nz, loop_dir, xbgn, xwid , rho):
  Nmx = nmp.max([Nx, Ny, Nz])

  rh = nmp.zeros(Nmx,'cdouble')
  rhstd = nmp.zeros(Nmx,'cdouble')
  nrm = nmp.zeros(Nmx, 'd')
  
  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz

  avg_shift = 0.0;
  for i in range( 0 , L1*L2 ):
    avg_shift += xbgn[i][1] - xbgn[i][0]

  avg_shift = int( avg_shift / L1 / L2 )


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      minim = xbgn[ind2][1] - avg_shift - 10
      if (minim < 0):
        minim = 0

      maxim = xbgn[ind2][1] + 10
      if (maxim >= LL):
        maxim = LL

      for j in range(minim, maxim):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rh[ j-minim ] += rho[ind]
 
        nrm[ j-minim ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rh[i] = rh[i] / nrm[i]


  for k in range(0 , L1):

    for i in range(0, L2):
 
      ind2 = k + i*L1
      
      minim = xbgn[ind2][1] - avg_shift - 7
      if (minim < 0):
        minim = 0
      
      maxim = xbgn[ind2][1] + 7
      if (maxim >= LL):
        maxim = LL


      for j in range(minim, maxim):

        ind = 0
 
        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        rhstd[ j-minim ] += (rh[ j-minim ] - rho[ind])**2
 
        nrm[ j-minim ] += 1.0


  for i in range(0, LL):
    if (nrm[i] > 0.0):
      rhstd[i] = (rhstd[i] / nrm[i]) ** 0.5


  return rh , rhstd

def find_3D_middle(x, y, z, loop_dir , thrhld, rho):
  Nx = nmp.shape(x)[0]
  Ny = nmp.shape(y)[0]
  Nz = nmp.shape(z)[0]

  dx = x[1] - x[0]
  
  Mnl = Nx * Ny

  xbgn = nmp.zeros((Mnl,2), 'int')
  xmid = nmp.zeros(Mnl , 'int' )

  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz


  for k in range(0 , L1):

    for i in range(0, L2):

      bgn_ind = k + i*L1

      for j in range(0, LL):


        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

        if ( xbgn[bgn_ind][0]==0 and  rho[ind] >= thrhld ):
          xbgn[ bgn_ind][0] = j


        if ( xbgn[ bgn_ind ][0] != 0 and rho[ind] < thrhld):
          xbgn[ bgn_ind ][1] = j
          xmid[ bgn_ind ] = (xbgn[bgn_ind][1] - xbgn[bgn_ind][0]) 
          break ;
 

        if (j==(LL-1)):
          print "interface not found!"
          exit(1)

  return xbgn , xmid



def find_3D_interface_indices(x, y, z, loop_dir , xcent, thrhld, rho):
  Nx = nmp.shape(x)[0]
  Ny = nmp.shape(y)[0]
  Nz = nmp.shape(z)[0]

  dx = x[1] - x[0]
  xci = int( xcent/dx )
  
  Mnl = Nx * Ny

  xbgn = nmp.zeros(Mnl, 'int')

  L1 = L2 = LL = 0;

  if (loop_dir == 0):
    LL = Nx
    L1 = Ny
    L2 = Nz
  elif (loop_dir == 1):
    L1 = Nx
    LL = Ny
    L2 = Nz
  else:
    L1 = Nx
    L2 = Ny
    LL = Nz


  for k in range(0 , L1):

    for i in range(0, L2):

      bgn_ind = k + i*L1

      for j in range(xci, LL):

        xbgn[i] = xci
 

        if (loop_dir==0):
          ind = j + (k + i*Ny)*Nx
        elif (loop_dir == 1):
          ind = k + (j + i*Ny)*Nx
        else:
          ind = k + (i + j*Ny)*Nx

 
        if (rho[ind] < thrhld):
          xbgn[ bgn_ind ] = j - 10 
          break ;
 
        if (j==(LL-1)):
          print "interface not found!"
          exit(1)

  return xbgn

def find_2D_interface_indices(x, y, xcent, thrhld, rho):
  Nx = nmp.shape(x)[0]
  Ny = nmp.shape(y)[0]

  dx = x[1] - x[0]
  xci = int( xcent/dx )
  xbgn = nmp.zeros(Ny, 'int')

  for i in range(0, Ny):
    xbgn[i] = xci

    for j in range(xci, Nx):
      ind = j + Nx*i

      if (rho[ind] < thrhld):
        xbgn[i] = j - 20
        break ;

      if (j==(Nx-1)):
        print "interface not found!"
        exit(1)

  return xbgn

def avg_dens_slice(Nx, Ny, xbgn, rho):

  rh = nmp.zeros(Nx,'cdouble')
  nrm = nmp.zeros(Ny, 'd')

  for i in range(0, Ny):

    for j in range(xbgn[i], Nx):

      ind = j + Ny*i

      rh[ j-xbgn[i] ] += rho[ind]

      nrm[ j-xbgn[i] ] += 1.0

  for i in range(0, Nx):
    if (nrm[i] > 0.0):
      rh[i] = rh[i] / nrm[i]

  return rh

def read_2d_density(Nx, Ny, name):

  M = Nx*Ny

  x = nmp.zeros(Nx,'d')
  y = nmp.zeros(Ny,'d')
  rho = nmp.zeros(M, 'cdouble')

  inp = open(name, 'r')

  line = inp.readline().split()
  for iy in range(0, Ny):
    for ix in range(0,Nx):
      id = ix + iy*Nx
       
      rho[id] = float(line[2]) + 1j*float(line[3])

      x[ix] = float(line[0])
      y[iy] = float(line[1])

      line = inp.readline().split()

    line = inp.readline().split()

  inp.close()

  return x, y, rho

def write_2d_density(name, x, y, dat):
  Nx = nmp.shape(x)[0]
  Ny = nmp.shape(y)[0]

  otp = open(name,'w')
  for iy in range(0, Ny):
    for ix in range(0, Nx):
      ind = ix + iy * Nx
      line = "%lf %lf %lf %lf\n" % (x[ix], y[iy], dat[ind].real, dat[ind].imag)
      otp.write(line)

    otp.write('\n')

  otp.close()

def read_real_2d_dat(Nx, Ny, name):
  M = Nx*Ny
  
  x = nmp.zeros(Nx,'double')
  y = nmp.zeros(Ny,'double')
  dat = nmp.zeros(M,'cdouble')

  inp = open(name,'r')
  ind = 0
  for iy in range(0,Ny):

    for ix in range(0, Nx):
      line = inp.readline().split()

      if (len(line) != 3):
        print "Line length error!"
        exit(1)

      x[ix] = float(line[0])
      y[iy] = float(line[1])

      dat[ind] = float(line[2]) 
      ind += 1
    line = inp.readline().split()
  inp.close()

  return x,y, dat

def read_2d_dat(Nx, Ny, name):
  M = Nx*Ny
  
  x = nmp.zeros(Nx,'double')
  y = nmp.zeros(Ny,'double')
  dat = nmp.zeros(M,'cdouble')

  inp = open(name,'r')
  ind = 0
  for iy in range(0,Ny):

    for ix in range(0, Nx):
      line = inp.readline().split()

      if (len(line) != 4):
        print "Line length error!"
        exit(1)

      x[ix] = float(line[0])
      y[iy] = float(line[1])

      dat[ind] = float(line[2]) + 1j*float(line[3])
      ind += 1
    line = inp.readline().split()
  inp.close()

  return x,y, dat

def read_real_3d_dat(Nx, Ny, Nz, name, rho_col=3, start_line=0):
  M = Nx*Ny*Nz
  
  x = nmp.zeros(Nx,'double')
  y = nmp.zeros(Ny,'double')
  z = nmp.zeros(Nz,'double')
  dat = nmp.zeros(M,'double')

  inp = open(name,'r')
  # skip headers
  for iline in range(start_line): inp.readline()

  ind = 0
  for iz in range(0, Nz):

    for iy in range(0,Ny):

      for ix in range(0, Nx):
        line = inp.readline().split()
        
        if (len(line) < rho_col):
          print "Line length error!"
          exit(1)

        x[ix] = float(line[0])
        y[iy] = float(line[1])
        z[iz] = float(line[2])

        dat[ind] = float(line[rho_col]) 
        ind += 1
  inp.close()

  return x,y,z, dat

def read_3d_dat(Nx, Ny, Nz, name):
  M = Nx*Ny*Nz
  
  x = nmp.zeros(Nx,'double')
  y = nmp.zeros(Ny,'double')
  z = nmp.zeros(Nz,'double')
  dat = nmp.zeros(M,'cdouble')

  inp = open(name,'r')
  ind = 0
  for iz in range(0, Nz):

    for iy in range(0,Ny):

      for ix in range(0, Nx):
        line = inp.readline().split()
        if (len(line) != 5):
          print "Line length error!"
          exit(1)

        x[ix] = float(line[0])
        y[iy] = float(line[1])
        z[iz] = float(line[2])

        dat[ind] = float(line[3]) + 1j*float(line[4])
        ind += 1
  inp.close()

  return x,y,z, dat

def read_orien_dat(Nx, Ny, Nz, name):
  M = Nx*Ny*Nz
  
  x = nmp.zeros(Nx,'double')
  y = nmp.zeros(Ny,'double')
  z = nmp.zeros(Nz,'double')
  datxx = nmp.zeros(M,'cdouble')
  datyy = nmp.zeros(M,'cdouble')
  datzz = nmp.zeros(M,'cdouble')
  datxy = nmp.zeros(M,'cdouble')
  datxz = nmp.zeros(M,'cdouble')
  datyz = nmp.zeros(M,'cdouble')

  inp = open(name,'r')
  ind = 0
  for iz in range(0, Nz):

    for iy in range(0,Ny):

      for ix in range(0, Nx):
        line = inp.readline().split()

        x[ix] = float(line[0])
        y[iy] = float(line[1])
        z[iz] = float(line[2])

        datxx[ind] = float(line[3])
        datyy[ind] = float(line[4])
        datzz[ind] = float(line[5])
        datxy[ind] = float(line[6])
        datxz[ind] = float(line[7])
        datyz[ind] = float(line[8])
        ind += 1
  inp.close()

  return x,y,z, datxx, datyy, datzz, datxy, datxz, datyz

def write_3d_dat(x, y, z, dat, name):
  Nx = nmp.shape(x)[0]
  Ny = nmp.shape(y)[0]
  Nz = nmp.shape(z)[0]

  M = Nx*Ny*Nz
  
  otp = open(name, 'w')
  ind = 0
  for iz in range(0, Nz):
    for iy in range(0, Ny):
      for ix in range(0, Nx):

        line = '%lf %lf %lf %e %e\n' % \
            (x[ix], y[iy], z[iz], dat[ind].real, dat[ind].imag)

        otp.write(line)
        ind += 1

  otp.close()

