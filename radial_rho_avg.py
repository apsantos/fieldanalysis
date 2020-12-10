#!/usr/bin/python
import math as m
import os, sys
import numpy as np

configs = ['100','200','300','400','500']

if ( len( sys.argv ) < 3 ):
    print "Usage: radial_avg_rho.py [dr_bin] [path] optional: [x posit] [y posit] [z posit]"

else:

    drbin = float( sys.argv[1] )
  
    path = sys.argv[2]
    
    for c in configs:
    
        nameA = path+c+'/avg_rhoga.dat'
        nameB = path+c+'/avg_rhogb.dat'
        
        print(nameA)
        
        inpA = open(nameA, 'r' )
        
        Dim = 3
        
        line = inpA.readline().split()
        
        ndat = 0
        while ( len( line ) > 0):
        
            ndat += 1
            line = inpA.readline().split()
        
    
        inpA.close()

        print ndat
        
        inpA = open( nameA , 'r' )
        inpB = open( nameB , 'r' )
        
        rhoA  = np.zeros( (ndat,Dim+1) , 'd' )
        rhoB  = np.zeros( (ndat,Dim+1) , 'd' )
        
        rhoAt  = np.zeros( (ndat,Dim+1) , 'd' )
        rhoBt  = np.zeros( (ndat,Dim+1) , 'd' )
        
        lineA = inpA.readline().split()
        lineB = inpB.readline().split()
    
        for i in range( 0 , ndat ):
        
            for j in range( 0 , Dim+1 ):
                rhoA[i,j] = float( lineA[j] )
                rhoB[i,j] = float( lineB[j] )
        
            lineA = inpA.readline().split()
            lineB = inpB.readline().split()
            
            if ( len( lineA ) == 0 ):
                lineA = inpA.readline().split()
                
            if ( len( lineB ) == 0 ):
                lineB = inpB.readline().split()
        
        inpA.close()
        inpB.close()

        Lmax = np.zeros( Dim , 'd' )
        for j in range( 0 , Dim ):
            Lmax[j] = np.max( rhoA[:,j] )
        
        nbins = int( np.max( Lmax )/drbin ) 
        
        print nbins, "bins allocated"
        
        norm = np.zeros( nbins, 'd' )
        radA = np.zeros( nbins , 'd' )
        radB = np.zeros( nbins , 'd' )
        
        radAt = np.zeros( nbins , 'd' )
        radBt = np.zeros( nbins , 'd' )
        
        cent = Lmax / 2.0 
        
        if ( len( sys.argv ) >= 5 ):
            cent[0] = float( sys.argv[3] )
            cent[1] = float( sys.argv[4] )
        
        if ( len( sys.argv ) == 6 ):
            cent[2] = float( sys.argv[5] )
        
        print "particle position: " , cent

        for i in range( 0 , ndat ):
            mdr2 = 0.0
            for j in range( 0 , Dim ):
                dr = rhoA[i,j] - cent[j]
                mdr2 += (rhoA[i,j] - cent[j]) ** 2
        
            mdr = m.sqrt( mdr2 )
        
            cur_bin = int( mdr / drbin )
        
            if ( cur_bin < nbins ):
                norm[ cur_bin ] += 1
                radA[ cur_bin ] += rhoA[i,Dim]
            
                radAt[ cur_bin ] += rhoA[i,Dim]
        
        outA = path+c+'/radial_rhoA.dat'
        otp = open( outA , 'w' )
        
        for i in range( 0 , nbins ):
            if ( norm[i] > 0.0 ):
                line = '%f %f\n' % ( (i+0.5) * drbin , radA[i]/norm[i] )
            else:
                line = '%f 0.0\n' % ( (i+0.5) * drbin )
        
            otp.write( line )
        
        otp.close()
        
        # B profiles
        for i in range( 0 , ndat ):
            mdr2 = 0.0
            for j in range( 0 , Dim ):
                dr = rhoB[i,j] - cent[j]
                mdr2 += (rhoB[i,j] - cent[j]) ** 2
        
            mdr = m.sqrt( mdr2 )
        
            cur_bin = int( mdr / drbin )
        
            if ( cur_bin < nbins ):
                norm[ cur_bin ] += 1
                radB[ cur_bin ] += rhoB[i,Dim]
            
                radBt[ cur_bin ] += rhoB[i,Dim]
        
        outB = path+c+'/radial_rhoB.dat'
        otp = open( outB , 'w' )
        
        for i in range( 0 , nbins ):
            if ( norm[i] > 0.0 ):
                line = '%f %f\n' % ( (i+0.5) * drbin , radB[i]/norm[i] )
            else:
                line = '%f 0.0\n' % ( (i+0.5) * drbin )
        
            otp.write( line )
        
        otp.close()
        
    outA = path+c+'/radial_rhoAt.dat'
    otp = open( outA, 'w')
    for i in range ( 0, nbins):
        if ( norm[i] > 0.0 ):
            line = '%f %f\n' % ( (i+0.5) * drbin , radAt[i]/(norm[i]*5) )
        else:
            line = '%f 0.0\n' % ( (i+0.5) * drbin )
        
        otp.write( line )
    otp.close()
    
    outB = path+c+'/radial_rhoBt.dat'
    otp = open( outB, 'w')
    for i in range ( 0, nbins):
        if ( norm[i] > 0.0 ):
            line = '%f %f\n' % ( (i+0.5) * drbin , radBt[i]/(norm[i]*5) )
        else:
            line = '%f 0.0\n' % ( (i+0.5) * drbin )
        
        otp.write( line )
    otp.close()
  
  
    
print('done')