#!/usr/bin/env python
"""np_properties

 Post processing particle positions of
 polymer nano-composites
 							                                    
 -Andrew P. Santos					                            

"""
import sys, argparse
from trajectory import Trajectory
import numpy as np
from numpy import linalg as LA
import scipy as sci
from scipy import optimize
import math as ma        
import time
import collections

compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class Analysis(object):
    """Callable Analyze class

    Run property calculation from input files using other classes

    """
    #def __init__(self):

    def addParser(self, parser):
        """
        Get relevant values from an argparse parser
        """
        self.start_frame = 1
        if (parser.parse_args().start_frame):
            self.start_frame = parser.parse_args().start_frame

        self.end_frame = 1
        if (parser.parse_args().end_frame):
            self.end_frame = parser.parse_args().end_frame

        self.ex_vol_calc = False
        if (parser.parse_args().ex_vol):
            self.ex_vol_calc = True
            self.ex_vol_n_conf = int(parser.parse_args().ex_vol[0])
            if (len(parser.parse_args().ex_vol) == 1):
                print 'Assuming user wants to use the saved configs'
                self.ex_vol_file = None

            elif (parser.parse_args().ex_vol[1] == 'saved'):
                self.ex_vol_file = None

            elif ('xyz' in parser.parse_args().ex_vol[1]):
                self.ex_vol_file = parser.parse_args().ex_vol[1]

            else:
                print 'CANNOT PERFORM EXCLUDED VOLUME CALCULATION'
                print 'Only xyz files or saved are options for '
                print 'Excluded Volume monomer configurations'
                self.ex_vol_calc = False

        if (parser.parse_args().ex_method):
            self.ex_method = parser.parse_args().ex_method
        else:
            self.ex_method = 'overlap'

        self.gyration_calc = parser.parse_args().gyration
        self.principal_calc = parser.parse_args().principal
        self.moment_calc = parser.parse_args().moment
        self.box = parser.parse_args().box_length

        self._traj_filename = parser.parse_args().traj_file[0]
        self._f_type = parser.parse_args().traj_type
        self._run_name = self._traj_filename[:len(self._traj_filename)-1-len(self._f_type)]

        self.np_type = parser.parse_args().np_type
        self.graft_type = parser.parse_args().graft_type
        self.matrix_type = parser.parse_args().matrix_type

        self.output_on = False
        if (parser.parse_args().output):
            self.output_on = True
            self.output_type = parser.parse_args().output
            if (self.output_type == 'frag'):
                if (not parser.parse_args().mcf):
                    print 'CANNOT CREATE Frag library'
                    print 'Requires lj sigma information, give mcf or ljparam file'
                    self.output_on = False
                elif (not parser.parse_args().box_length):
                    print 'CANNOT CREATE Frag library'
                    print 'Requires box size to get the COM correction right'
                    self.output_on = False
            elif (self.output_type == 'xml'):
                if (parser.parse_args().xml_d):
                    self.xml_d = parser.parse_args().xml_d
                else:
                    print 'CANNOT WRITE XML file'
                    print 'Requires d_scale'
                    self.xml_d = 0.0

                if (parser.parse_args().xml_e):
                    self.xml_e = parser.parse_args().xml_e
                else:
                    self.xml_e = 1.0

                if (parser.parse_args().xml_m):
                    self.xml_m = parser.parse_args().xml_m
                else:
                    self.xml_m = 1.0

                self.special_pairs = []
                if (parser.parse_args().special_pair):
                    n_special_pair = int( len(parser.parse_args().special_pair) )
                    for i in range(n_special_pair):
                        self.special_pairs.append( parser.parse_args().special_pair[i] )
                    

        self.density_calc = False
        if (parser.parse_args().density):
            self.density_calc = True
            self.density_bin_width = parser.parse_args().density

        self.rdf_calc = False
        self.rdf_type_calc = False
        self.rdf_inNP_calc = False
        if (parser.parse_args().rdf or parser.parse_args().rdf_inNP):
            self.rdf_calc = True
            self.rdf_bin_width = [0, 0, 0]
            if (parser.parse_args().rdf):
                self.rdf_type_calc = True
                self.rdf_bin_width[0] = float(parser.parse_args().rdf[0])
                n_rdf_types = int( (len(parser.parse_args().rdf) - 1) / 2.0)
                self.rdf_types = np.chararray( (n_rdf_types, 2), itemsize=3)
                i = 1
                itype = 0
                while itype < n_rdf_types:
                    self.rdf_types[itype,0] = parser.parse_args().rdf[i]
                    self.rdf_types[itype,1] = parser.parse_args().rdf[i+1]
                    i += 2
                    itype += 1

            elif (parser.parse_args().rdf_inNP):
                self.rdf_inNP_calc = True
                self.rdf_bin_width[1] = float(parser.parse_args().rdf_inNP[0])
                n_rdf_types = int( (len(parser.parse_args().rdf_inNP) - 1) / 2.0)
                self.rdf_types = np.chararray( (n_rdf_types, 2), itemsize=3)
                i = 1
                itype = 0
                while itype < n_rdf_types:
                    self.rdf_types[itype,0] = parser.parse_args().rdf_inNP[i]
                    self.rdf_types[itype,1] = parser.parse_args().rdf_inNP[i+1]
                    i += 2
                    itype += 1

            else:
                if not self.rdf_type_calc:
                    self.rdf_calc = False
                self.rdf_inNP_calc = False

    def initialize(self):
        """
        Add analysis values based on a string
        """
        if (self.density_calc):
            self.density = Density(self.density_bin_width, self.param_name)

        if (self.moment_calc or self.gyration_calc):
            self.RG = RadiiGyration()
            self.RG.moment_calc = self.moment_calc
            self.RG.gyration_calc = self.gyration_calc
            self.RG.principal_calc = self.principal_calc

        if (self.rdf_calc):
            maxM = 0
            self.rdf = rdf(self.rdf_bin_width, maxM)
            self.rdf.type_calc = self.rdf_type_calc
            self.rdf.inNP_calc = self.rdf_inNP_calc
            if self.rdf_type_calc or self.rdf_inNP_calc:
                self.rdf.setTypes(self.rdf_types)

            self.rdf.setBox(self.box)
            self.rdf.setBinProperties()

    def analyze(self, trajectory):
        """
        Run the analysis of the analyzers added
        """
        trajectory.open()
        self.initialize()
        #n_frames = trajectory.getNframes()
        n_frames = trajectory._n_frames
        if n_frames < self.start_frame:
            print 'There are only', n_frames, 'frames in the trajectory'
            print self.start_frame, 'frames were given'
            return 1
        elif self.end_frame == 0:
            self.end_frame = n_frames
        elif self.end_frame < self.start_frame:
            print ('The starting frame number (', self.start_frame, 
                  ') must be smaller than the ending frame (', self.end_frame, ')')
            return 1

        print 'Total frames in trajectory:\t', n_frames
        print 'Number of frames processing:\t', (self.end_frame - self.start_frame + 1)
        sys.stdout.write('reading frame: ')
        sys.stdout.flush()
        for iframe in range(1,n_frames+1):
            frame = trajectory.read(iframe < self.start_frame)
            # Decide whether to use frame or not
            if iframe < self.start_frame:
                continue

            elif iframe > self.end_frame:
                break

            elif (iframe == n_frames):
                sys.stdout.write(str(n_frames)+'\n')
                sys.stdout.flush()

            elif (iframe % 100 == 0 or iframe == self.start_frame):
                sys.stdout.write(str(iframe)+', ')
                sys.stdout.flush()

            if (self.output_on):
                trajectory.output_on = True
                if iframe == self.start_frame:
                    # set the box dimensions
                    trajectory.setBox(self.box)
                    self.setTrajProperties(frame, trajectory)

            if self.rdf_inNP_calc or self.density_calc:
                NPmolids = self.find_NPmolids(frame)
                for iNP in NPmolids:
                    if self.rdf_type_calc:
                        if self.rdf_inNP_calc:
                            self.rdf.addNP(iNP)
                    if self.density_calc:
                        self.density.calculate(iNP)

            if self.rdf_type_calc:
                self.rdf.setFrame(frame)
                self.rdf.calculateType()
                if self.rdf_inNP_calc:
                    self.rdf.calculateInNP()


            if (self.output_on):
                if (self.savethemicelles):
                    # Don't print empty clusters
                    if (frame['natoms'] == 0):
                        continue

                trajectory.write(frame)

                if (self.output_type == 'chk'):
                    sys.stdout.write('Output checkpoint file')
                    break

        sys.stdout.write('\nCalculating the')
        if (self.density_calc):
            sys.stdout.write(', density profile')

        if (self.rdf_calc):
            sys.stdout.write(', RDFs')

        sys.stdout.write('\n')
        sys.stdout.flush()

        trajectory.close()

    def write(self):
        """
        Write the analysis results
        """
        if self.density_calc:
            self.density.write(self._run_name)
        if self.gyration_calc or self.moment_calc:
            self.RG.write(self._run_name)
        if self.rdf_calc:
            if self.rdf_type_calc:
                self.rdf.writeType(self._run_name)

    def findNPmolids(self, frame):
        molids = []
        for iatom in range(frame['natoms']):
            if frame['aname'][iatom] == self.NPname:
                molids.append( frame['amol'][iatom] )

        return molids 

class rdf(object):
    """
    Calculates and writes the radial distribution function
    """
    def __init__(self, binsize, maxclus=20):
        self.binsize_type = binsize[0]
        self.binsize_inNP = binsize[1]
        if self.binsize_type <= 0 and self.binsize_inNP <= 0: 
            print 'cannot calculate rdf binsizes need to be > 0'
        else:
            self.binsize = 0.0
            for ibinsize in binsize:
                self.binsize = max(self.binsize, ibinsize)
        
        self.nNP = 0
        self.bin_cut = 0
        self.type_calc = False
        self.inNP_calc = False

    def setTypes(self, types):
        self.types = types
        self.n_type_groups = len(self.types[:,0])
        self.n_total_group = np.zeros( (self.n_type_groups), dtype = np.int)
        self.nconf_inNP = np.zeros( (self.n_type_groups), dtype = np.int)
        self.nconf_type = np.zeros( (self.n_type_groups), dtype = np.int)

    def setBox(self, box_length):
        # assume the lower bound is 0
        if (len(box_length) == 3 or len(box_length) == 2):
            if (len(box_length) == 3):
                self.dim = 3
            elif (len(box_length) == 2):
                self.dim = 2
            self.box = np.zeros( (self.dim, 2) )
            self.box_length = np.array(box_length)
            self.box_length_half = np.array(box_length) / 2.0
            self.box_mid = np.array(box_length) / 2.0
            for idim in range(self.dim):
                self.box[idim, 1] = box_length[idim]

        # user gave both upper and lower bound
        elif (len(box_length) == 6 or len(box_length) == 4):
            if (len(box_length) == 6):
                self.dim = 3
            elif (len(box_length) == 4):
                self.dim = 2
            self.box = np.zeros( (self.dim, 2) )
            self.box_length = np.zeros((self.dim))
            self.box_length_half = np.zeros((self.dim))
            self.box_mid = np.zeros((self.dim))
            for idim in range(self.dim):
                for ilim in range(2):
                    self.box[idim, ilim] = box_length[idim*2+ilim]

                self.box_length[idim] = self.box[idim, 1] - self.box[idim, 0]
                self.box_length_half[idim] = (self.box[idim, 1] - self.box[idim, 0]) / 2.0
                self.box_mid[idim] = (self.box[idim, 1] + self.box[idim, 0]) / 2.0
                if (self.box[idim, 1] < self.box[idim, 0]):
                    raise ValueError('Box ill defined, '
                                 'lower bound must be lower than upper bound')
        else:
            raise TypeError('Box ill, or not defined, '
                           'not enough entrants (3 or 6 for 3D, 2 or 4 for 2D')

    def setFrame(self, frame):
        self.natom = frame['natoms']
        self.aname = frame['aname']
        self.pos = frame['pos']
        self.amol = frame['amol']

        self.setTypeSizes()

    def addNP(self, imolNP):
        #cluster class
        self.imolNP = imolNP
        if self.inNP_calc:
            self.setTypeSizes()

    def setTypeSizes(self):
        for igroup in range(self.n_type_groups):
            for iatom in range(self.natom):
                if (self.aname[iatom] == self.types[igroup,0]):
                    self.n_total_group[igroup] += 1
                if (self.aname[iatom] == self.types[igroup,1]):
                    self.n_total_group[igroup] += 1
    
    def setBinProperties(self):
        self.bin_cut = min(self.box_length) / 2.0
        self.bin_cut2 = self.bin_cut**2.0

        self.nbins = int(ma.ceil(self.bin_cut / self.binsize))

        if self.type_calc:
            self.nbins_type = int(ma.ceil(self.bin_cut / self.binsize_type))
            self.gr_type = np.zeros( (self.nbins_type+1, len(self.types[:,0])), dtype = np.float)

        if self.inNP_calc:
            self.nbins_inNP = int(ma.ceil(self.bin_cut / self.binsize_inNP))
            self.gr_inNP = np.zeros( (self.nbins_inNP+1, len(self.types[:,0])), dtype = np.float)

        self.vol = 1.0
        for i in range(len(self.box_length)):
            self.vol *= self.box_length[i]

        self.id_prefactor = (4.0/3.0) * ma.pi / self.vol

    def writeType(self, runname):
        for igroup in range(self.n_type_groups):
            ofile = open('./%s.%s-%s.rdf' % (runname, self.types[igroup,0], self.types[igroup,1]), 'w')
            ofile.write( "r[A]   rdf\n")
            ofile.write( "%f %f\n" % (0.0, 0.0) )
                
            factor = self.id_prefactor * self.nconf_type[igroup]
            #if (self.types[igroup,0] != self.types[igroup,1]):
            #    factor *= 2.0
            for ibin in range(self.nbins_type):
                r = self.binsize_type * (ibin + 0.5)
                V_bin = self.binsize_type**3.0 * ((ibin + 1)**3.0 - ibin**3.0)
                self.gr_type[ibin, igroup] /= float(factor * V_bin)
                ofile.write( "%f %f\n" % (r, self.gr_type[ibin, igroup]) )
            ofile.close()

    def writeInNP(self, runname):
        for igroup in range(self.n_type_groups):
            ofile = open('./%s.%s-%s.inNP.rdf' % (runname, self.types[igroup,0], self.types[igroup,1]), 'w')
            ofile.write( "r[A]   rdf\n")
            if self.ncluster == 0:
                ofile.write( "no Nanoparticles found" )
                return
            
            ofile.write( "%f %f\n" % (0.0, 0.0) )
                
            factor = self.id_prefactor * self.nconf_inNP[igroup]
            #if (self.types[igroup,0] != self.types[igroup,1]):
            #    factor *= 2.0
            for ibin in range(self.nbins_inNP):
                r = self.binsize_inNP * (ibin + 0.5)
                V_bin = self.binsize_inNP**3.0 * ((ibin + 1)**3.0 - ibin**3.0)
                self.gr_inNP[ibin, igroup] /= float(factor * V_bin)
                ofile.write( "%f %f\n" % (r, self.gr_inNP[ibin, igroup]) )
            ofile.close()

    def calculateType(self):
        group_check = np.zeros( (self.n_type_groups), dtype = np.int)
        for igroup in range(self.n_type_groups):
            # make a list of all atoms in the cluster
            i_atoms = []
            j_atoms = []
            for iatom in range(self.natom):
                if (self.aname[iatom].strip() == self.types[igroup,0]):
                    i_atoms.append(iatom)
                if (self.aname[iatom].strip() == self.types[igroup,1]):
                #if (self.aname[iatom] == self.types[igroup,1]):
                    j_atoms.append(iatom)

            n_i = len(i_atoms)
            n_j = len(j_atoms)
            if (n_i == 0 or n_j == 0): 
                continue
            self.nconf_type[igroup] += 1

            same_type = (self.types[igroup,0] == self.types[igroup,1])

            t_gr = self.calcRDFtypes(i_atoms, j_atoms, self.nbins_type, self.binsize_type, same_type)
            for ibin in range(self.nbins_type):
                self.gr_type[ibin, igroup] += t_gr[ibin] / n_i / n_j 

    def calculateInNP(self):
        group_check = np.zeros( (self.n_type_groups), dtype = np.int)
        for iclus in self.l_cluster:
            iM = self.N[ iclus ]
            if (iM in self.clus_list):

                self.ncluster += 1

                for igroup in range(self.n_type_groups):
                    # make a list of all atoms in the cluster
                    i_clus_atoms = []
                    j_clus_atoms = []
                    for iatom in range(self.natom):
                        # if the amphiphile is in the cluster
                        if (self.clabel[self.amol[iatom]] == iclus):
                            if (self.aname[iatom].strip() == self.types[igroup,0]):
                                i_clus_atoms.append(iatom)
                            if (self.aname[iatom].strip() == self.types[igroup,1]):
                                j_clus_atoms.append(iatom)
        
                    n_i = len(i_clus_atoms)
                    n_j = len(j_clus_atoms)
                    if (n_i == 0 or n_j == 0): 
                        continue
                    #elif group_check[igroup] == 0:
                    #    self.nconf_inNP[igroup] += 1
                    #    group_check[igroup] = 1
                    self.nconf_inNP[igroup] += 1

                    if (self.types[igroup,0] == self.types[igroup,1]):
                        same_type = True
                    else:
                        same_type = False

                    t_gr = self.calcRDFtypes(i_clus_atoms, j_clus_atoms, self.nbins_inNP, self.binsize_inNP, same_type)
                    for ibin in range(self.nbins_inNP):
                        self.gr_inNP[ibin, igroup] += t_gr[ibin] / n_i / n_j 
                        #self.gr_inNP[ibin, igroup] += t_gr[ibin] / ma.sqrt( n_i * n_j )

    def calcRDFtypes(self, i_list, j_list, nbin, binsize, same_type=True):
        t_gr = np.zeros( (nbin+1), dtype = np.float)

        if same_type:
            i_tmp = i_list[:len(i_list)-1]
            i = 0
        else:
            i_tmp = i_list

        for iatom in i_tmp:
            if (same_type):
                j_tmp = i_list[i+1:]
                i += 1
            else:
                j_tmp = j_list

            for jatom in j_tmp:
                if (self.amol[iatom] == self.amol[jatom]): continue

                r2 = 0
                for idim in range(self.dim):
                    # calculate the positions w.r.t. the COM
                    d = self.pos[iatom, idim] - self.pos[jatom, idim]
                    # correct for if the position is cut by the boundary
                    dist = d - (self.box_length[idim] * int( round(d / self.box_length[idim]) ))
                    r2 += dist**2.0

                if (r2 < self.bin_cut2):
                    ibin = int(r2**0.5 / binsize)
                    #if ibin < 14:
                    #    print ibin, r2
                    if (same_type):
                        t_gr[ibin] += 2
                    else:
                        t_gr[ibin] += 1
            
        return t_gr

class RadiiGyration(object):
    """
    Calculates and writes the radii of gyration
    """
    def __init__(self, Mmax=20):
        self.max_NP = Mmax
        self.Rg = np.zeros( (self.max_NP+1, 4), dtype = np.float)
        self.moment = np.zeros( (self.max_NP+1, 4), dtype = np.float)
        self.nNP = np.zeros( (self.max_NP+1), dtype = np.int)
        self.moment_calc = False
        self.gyration_calc = False

    def addNP(self, imolNP):
        self.imolNP = imolNP
        return

    def write(self, runname):
        if (sum(self.nNP) == 0):
            print 'No NPters found'
            return

        self.normalize()

        if self.gyration_calc:
            gy_ofile = open('./' + runname + '.radii', 'w')
            gy_ofile.write( "M principal lower last\n")

        if self.moment_calc:
            mom_ofile = open('./' + runname + '.moment', 'w')
            mom_ofile.write( "M total principal lower last\n")

        for iM in range(1, self.max_NP):
            if (self.nNP[iM] > 0):
                if self.gyration_calc:
                    gy_ofile.write( "%d %f %f %f %f\n" % (iM, self.Rg[iM,0], self.Rg[iM,1], self.Rg[iM,2], self.Rg[iM,3]) )
                if self.moment_calc:
                    mom_ofile.write( "%d %f %f %f %f\n" % (iM, self.moment[iM,0], self.moment[iM,1], self.moment[iM,2], self.moment[iM,3]) )

        if self.gyration_calc:
            gy_ofile.close()
    
            gy_ofile = open('./' + runname + 'ave.radii', 'w')
            gy_ofile.write( "principal lower last\n")
            gy_ofile.write( "%f %f %f %f\n" % (self.Rg_ave[0], self.Rg_ave[1], self.Rg_ave[2], self.Rg_ave[3]) )
            gy_ofile.close()

        if self.moment_calc:
            mom_ofile.close()
    
            mom_ofile = open('./' + runname + 'ave.moment', 'w')
            mom_ofile.write( "total principal lower last\n")
            mom_ofile.write( "%f %f %f %f\n" % (self.I_ave[0], self.I_ave[1], self.I_ave[2], self.I_ave[3]) )
            mom_ofile.close()
        
    def normalize(self):
        if (sum(self.ncluster) == 0):
            print 'No clusters found'
            return

        if self.gyration_calc:
            R_ave = np.zeros( (4), dtype = np.float)

        if self.moment_calc:
            moment_ave= np.zeros( (4), dtype = np.float)

        nclus_averaged = 0
        for iM in range(1, self.max_clus):
            if (self.ncluster[iM] > 0):
                if self.gyration_calc:
                    self.Rg[iM,:] = np.divide(self.Rg[iM], float(self.ncluster[iM]))
                    R_ave[0] += self.Rg[iM, 0]
                    R_ave[1] += self.Rg[iM, 1]
                    R_ave[2] += self.Rg[iM, 2]
                    R_ave[3] += self.Rg[iM, 3]
    
                if self.moment_calc:
                    self.moment[iM, :] = np.divide(self.moment[iM], float(self.ncluster[iM]))
                    moment_ave[0] += self.moment[iM, 0]
                    moment_ave[1] += self.moment[iM, 1]
                    moment_ave[2] += self.moment[iM, 2]
                    moment_ave[3] += self.moment[iM, 3]

                nclus_averaged += 1

        if self.gyration_calc:
            self.Rg_ave = np.divide(R_ave, float(nclus_averaged))

        if self.moment_calc:
            self.I_ave = np.divide(moment_ave, float(nclus_averaged))

    def calculate(self, M):
        this_clus = -1
        for iclus in self.clus.l_cluster:
            if (M == self.clus.N[iclus] and iclus not in self.lclusters):
                this_clus = iclus
                self.lclusters.append(iclus)
                break

        if (this_clus < 0): return
        # positions corrected for the COM
        com_clus, com_clus_pos = self.clus.comClusterPos(iclus)

        n_atom_cluster = len(com_clus_pos[:,0])
        if (n_atom_cluster < 2): return

        self.ncluster[M] += 1

        I_ij = np.zeros( (3, 3), dtype = np.float)
        if self.gyration_calc:
            x_ij = np.zeros( (3, 3), dtype = np.float)

        # loop over amphiphiles/chains/molecules
        for iatom in range( n_atom_cluster ):
            # calculate the moment of gyration tensor
            xx = com_clus_pos[iatom,0]**2
            yy = com_clus_pos[iatom,1]**2
            zz = com_clus_pos[iatom,2]**2
            I_ij[0,0] += ( yy + zz )
            I_ij[1,1] += ( xx + zz )
            I_ij[2,2] += ( xx + yy )
            I_ij[0,1] += -( com_clus_pos[iatom,0] * com_clus_pos[iatom,1] )
            I_ij[1,2] += -( com_clus_pos[iatom,1] * com_clus_pos[iatom,2] )
            I_ij[0,2] += -( com_clus_pos[iatom,0] * com_clus_pos[iatom,2] )
            if self.gyration_calc:
                x_ij[0,0] += xx
                x_ij[1,1] += yy
                x_ij[2,2] += zz

        I_ij[1,0] = I_ij[0,1]
        I_ij[2,1] = I_ij[1,2]
        I_ij[2,0] = I_ij[0,2]
        I_eigval = LA.eigvals(I_ij)
        # total moment
        self.moment[M, 0] += ma.sqrt(I_eigval[0]**2 + I_eigval[1]**2 + I_eigval[2]**2)
        # principle moments
        self.moment[M, 1] += np.amax(I_eigval)
        self.moment[M, 2] += np.median(I_eigval)
        self.moment[M, 3] += np.amin(I_eigval)
        if self.gyration_calc:
            self.Rg[M, 0] += ma.sqrt( (x_ij[0,0] + x_ij[1,1] + x_ij[2,2]) / n_atom_cluster )
            if self.principal_calc:
                x_ij[0,1] = -I_ij[0,1]
                x_ij[1,0] = -I_ij[1,0]
                x_ij[1,2] = -I_ij[1,2]
                x_ij[2,1] = -I_ij[2,1]
                x_ij[0,2] = -I_ij[0,2]
                x_ij[2,0] = -I_ij[2,0]
                x_eigval = LA.eigvals(x_ij)
                self.Rg[M, 1] += ma.sqrt( np.amax(x_eigval) / n_atom_cluster )
                self.Rg[M, 2] += ma.sqrt( np.median(x_eigval) / n_atom_cluster )
                self.Rg[M, 3] += ma.sqrt( np.amin(x_eigval) / n_atom_cluster )
            else:
                self.Rg[M, 1] += ma.sqrt( np.amax(x_ij) / n_atom_cluster )
                self.Rg[M, 2] += ma.sqrt( np.median(x_ij) / n_atom_cluster )
                self.Rg[M, 3] += ma.sqrt( np.amin(x_ij) / n_atom_cluster )
        
class Density(object):
    """
    Calculates and writes Density profile
    """
    def __init__(self, binsize, type_names, nbins=1000):
        self.binsize = binsize
        self.lclusters = []
        self.ncluster = 0
        self.bin_cut = 0
        self.ntypes = 0
        self.types = []
        self.nbins = nbins
        for ispecies in range(len(type_names)):
        #for ispecies in range(len(type_names[:][0])):
            for iname in type_names[ispecies][:]:
                if iname not in self.types:
                    self.types.append(iname)
                    self.ntypes += 1

        self.density = np.zeros( (self.ntypes, self.nbins), dtype = np.float)
        self.max_dist = 0

    def addCluster(self, cluster):
        #cluster class
        self.clus = cluster

        self.lclusters = []

    def write(self, runname):
        if (self.ncluster == 0):
            print 'No clusters found'
            return

        dens_ofile = open('./' + runname + '.rho', 'w')
        dens_ofile.write( "# r " )
        for iname in self.types:
            dens_ofile.write( "%s " % iname)
        dens_ofile.write( "\n" )
        factor = (4.0 / 3.0) * ma.pi
        for ibin in range(0, self.nbins):
            r = self.binsize * (ibin + 0.5)
            if (r > self.max_dist):break

            dens_ofile.write( "%f " % r)

            V_bin = 1#factor * self.binsize**3.0 * ((ibin + 1)**3.0 - ibin**3.0)
            for itype in range(self.ntypes):
                dens_ofile.write( "%f " % (self.density[itype, ibin] / self.ncluster / V_bin) )

            dens_ofile.write( "\n" )

    def calculate(self, M):
        this_clus = -1
        for iclus in self.clus.l_cluster:
            if (M == self.clus.N[iclus] and iclus not in self.lclusters):
                this_clus = iclus
                self.lclusters.append(iclus)
                break

        if (this_clus < 0): return
        # positions corrected for the COM
        com_clus, com_clus_pos, com_clus_name = self.clus.comClusterPosName(iclus)

        n_atom_cluster = len(com_clus_pos[:,0])
        if (n_atom_cluster < 2): return

        self.ncluster += 1
        # loop over amphiphiles/chains/molecules
        for iatom in range( n_atom_cluster ):
            # calculate the moment of gyration tensor
            distance = (com_clus_pos[iatom,0]**2 + com_clus_pos[iatom,1]**2 + com_clus_pos[iatom,2]**2)**0.5
            itype = 0
            for iname in self.types:
                if com_clus_name[iatom] == iname: break
                itype += 1

            if (distance > self.max_dist ) :
                self.max_dist = distance

            ibin = int(distance / self.binsize)
            self.density[itype, ibin] += 1

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Analyze configurations from simulations '
                                                 'of surfactant, chain molecules '
                                                 'which are likely aggregated')
    parser.add_argument('-t', "--traj_file", type=str, nargs='+',
                   help='Configuration/Trajectory/movie file: *xyz '
                        'Right now only 1 file capability, ' 
                        'pdb and towhee_movie are under development')
    parser.add_argument("--traj_type", type=str, choices=['xyz', 'gro', 'pdb'],
                   help='trajectory filetype')
    parser.add_argument("--start_frame", type=int,
                   help='Set the starting frame number.')
    parser.add_argument("--end_frame", type=int,
                   help='Set the ending frame number.')
    parser.add_argument('-d', "--density", type=float,
                   help='Calculate the density profile of each atom types in micelles, give the bin width')
    parser.add_argument('-I', "--moment", action="store_true",
                   help='Calculate the moments of inertia of micelles')
    parser.add_argument('-r', "--gyration", action="store_true",
                   help='Calculate the radius of gyration of micelles')
    parser.add_argument("--principal", action="store_true",
                   help='Calculate the radius of gyration of micelles around the principal axes')
    parser.add_argument('-e', "--ex_vol", type=str, nargs='+',
                   help='Number of monomer configurations for '
                        'excluded volume calculation, and method for monomer configs.'
                        'eg: -e energy 25 monomers.xyz')
    parser.add_argument("--ex_method", type=str, choices=['energy', 'overlap'],
                        help='method for excluded volume calculation either by the energy'
                             'or the overlap of sigmas. eg: --ex_method energy')
    parser.add_argument('-b', "--box_length", type=float, nargs='+',
                   help='Simulation box length (orthorombic is the only option)'
                        '0-box_lenth can give only 3')
    parser.add_argument("--rdf", type=str, nargs='+',
                   help='Radial distribution function, bin width and the atom names'
                        'eg: --rdf 0.3 IO HG')
    parser.add_argument("--rdf_inNP", type=str, nargs='+',
                   help='Radial distribution function within a nanoparticle bin width and the atom names'
                        'eg: --rdf_inNP 0.3 IO HG')
    parser.add_argument("--np_type", type=str, default='O',
                   help='charachter name for the NP type')
    parser.add_argument("--graft_type", type=str, default='H',
                   help='charachter name for the graft type')
    parser.add_argument("--matrix_type", type=str, default='He',
                   help='charachter name for the matrix type')
    parser.add_argument("-o", "--output", metavar='outputFile', type=str, choices=['chk', 'frag', 'xyz', 'gro', 'g96', 'xml', 'lmp'],
                   help='Output a file in a specific format (frag requires mcf to zero at COM):'
                        '-o xyz')

    # Initialize the Analysis class
    analyze = Analysis()

    # Tell the class everything specified in files and command line
    err = analyze.addParser(parser)

    if err != None:
        return

    # Initialize the Trajectory
    traj = Trajectory()
    traj.addFilename(parser.parse_args().traj_file[0], parser.parse_args().output)

    analyze_error = analyze.analyze(traj)

    if not analyze_error:
        analyze.write()

if __name__ == '__main__':
    sys.exit(main())

