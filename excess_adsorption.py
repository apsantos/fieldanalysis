#!/usr/bin/python
import sys, argparse
import numpy as np
from scipy import integrate

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit vapor, liquid density to get critical point and curve')
    parser.add_argument("--profile_file", type=str, required=True, help='input density profile file')
    parser.add_argument('-l', "--start_line", type=int, default=1, help='starting line of density profile file')
    parser.add_argument("--r_col", type=int, default=0, help='radius data collumn')
    parser.add_argument("--r_min", type=float, default=0.0, help='radius min to integrate')
    parser.add_argument("--r_max", type=float, default=1000000.0, help='radius min to integrate')
    parser.add_argument("--rho_col", type=int, default=1, help='density data collumn')
    parser.add_argument("--rho_inf", type=float, default=1.0, help='the far-field density')
    parser.add_argument("--gibbssurface", action='store_true', default=False, help='Find the gibbs dividing surface and use as the low r limit')
    c_r = parser.parse_args().r_col
    c_rho = parser.parse_args().rho_col
    rho_inf = parser.parse_args().rho_inf
    r_min = parser.parse_args().r_min
    r_max = parser.parse_args().r_max
    ndat = file_len(parser.parse_args().profile_file) - parser.parse_args().start_line

    # read in temeprature and vapor/liquid densities
    inp = open(parser.parse_args().profile_file, 'r') 

    j = 0
    r = []
    rho = []
    deltarho = []
    for line in inp:
        data = line.strip().split()
        if j < parser.parse_args().start_line:
            j += 1
            continue
        tr = float(data[c_r])
        if r_min < tr < r_max:
            r.append( tr )
            rho.append( float(data[c_rho]) )
            deltarho.append( float(data[c_rho])-rho_inf )

    r = np.array(r)
    deltarho = np.array(deltarho)
    if parser.parse_args().gibbssurface:
        otp = open(parser.parse_args().profile_file+'.derivative', 'w') 
        otp.write("# r rho smooth_rho drhodr integrand\n" )
        # get derivative of rho
        smooth_rho = smooth(rho,3)
        dr = r[1]-r[0]
        #drhodr = np.gradient(smooth_rho, 4*(r[1]-r[0]))
        drhodr = np.gradient(rho, dr)
        cutoff=0
        integrand = np.zeros((len(drhodr)-cutoff))
        for i in range(len(integrand)):
            V_bin = 3.*dr*r[i]**2. + dr**3./4.
            integrand[i] = rho[i]*V_bin
            #integrand[i] = drhodr[i]*V_bin
            #integrand[i] = drhodr[i]*r[i]**3.0
            otp.write("%f %f %f %f %f\n" % ( r[i], rho[i], smooth_rho[i], drhodr[i], integrand[i] ) )

        # liquid droplet tpye
        #re = (1/(0.0-rho_inf) * integrate.simps(integrand,x=r) )**(1./3.)
        # assuming it is away from a surface
        #re = ((1.0/rho_inf * integrate.simps(integrand,x=r[:len(r)-cutoff])))**(1./3.)
        #re = ( (1.0/rho_inf * integrate.simps(integrand,x=r[:len(r)-cutoff])))**(1./3.)
        #N = integrate.simps(integrand,x=r[:len(r)-cutoff])
        N = np.sum(integrand)
        Ve = max(0.0,r_max**3.0 - N/rho_inf)
        re = (Ve)**(1./3.)
        #print r_max, N
        for ire in range(len(r)):
            if r[ire] > re:
                break
        deltarho = deltarho[ire:]
        r = r[ire:]
        
        
    integral = integrate.simps(deltarho,x=r)
    
    if parser.parse_args().gibbssurface:
        print integral, re
    else:
        print integral

if __name__ == '__main__':
    sys.exit(main())

