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

        self.start_frame = parser.parse_args().start_frame
        self.end_frame = parser.parse_args().end_frame

        self.box = parser.parse_args().box_length
        self.box_length = parser.parse_args().box_length
        self.box_length_half = [0,0,0]
        for idim in range(3):
            self.box_length_half[idim] = self.box[idim]/2.0

        self._traj_filename = parser.parse_args().traj_file[0]
        self._f_type = parser.parse_args().traj_type

        self.center_type = parser.parse_args().center_type
        self.center = parser.parse_args().center
        self.np_aname = parser.parse_args().np_name
        self.graft_aname = parser.parse_args().graft_name
        self.anchor_aname = parser.parse_args().anchor_name
        self.matrix_aname = parser.parse_args().matrix_name

        self.gnp_mol_name = parser.parse_args().gnp_name
        self.solution_anames = [self.matrix_aname]

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
                    
    def analyze(self, trajectory):
        """
        Run the analysis of the analyzers added
        """
        trajectory.open()
        #self.initialize()
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
        for iframe in range(0,n_frames):
        #for iframe in range(1,n_frames+1):
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
                self.dim = 3
                if iframe == self.start_frame:
                    # set the box dimensions
                    trajectory.setBox(self.box)
                    #self.setTrajProperties(frame, trajectory)
                self.pos = frame['pos']
                self.atomids, self.molids = self.findCenterIds(frame)
                self.write(frame)
                frame['pos'] = self.shiftpos
                frame['amol'] = self.mols
                trajectory.write(frame)

        sys.stdout.write('\n')
        sys.stdout.flush()

        trajectory.close()

    def calcDist(self, ix, jx, iy=0.0, jy=0.0,iz=0.0,jz=0.0):
        dist = 0.0
        d = [(ix - jx), (iy - jy), (iz - jz)]
        for idim in [0,1,2]:
            #dist = d - (self.box_length[idim] * int( round(d / self.box_length[idim]) ))
            #elif (d[idim] >= self.box[idim, 1]):
            if   (d[idim] < -self.box_length_half[idim]):
                d[idim] += self.box_length[idim]
            elif   (d[idim] >= self.box_length_half[idim]):
                d[idim] -= self.box_length[idim]
            dist += d[idim]**2.0
        return dist


    def write(self,frame):
        """
        Write the shifted positions
        """
        # get the position to shift as function of molid
        self.shiftpos = np.zeros((len(self.pos[:,0]),3))
        self.mols = np.zeros((len(self.pos[:,0])))
        self.reference = np.zeros((3))
        nmols = len(self.molids)
        for imol in range(nmols):
            if self.center == 'graft':
                # ADD loop for NPS
                self.graftids, self.anchorids, self.npids = self.findGrafts(frame, imol)
                if self.reference[0] == 0:
                    self.R2 = self.calcDist(self.pos[self.anchorids[0],0], self.pos[self.npids[0],0],
                                            self.pos[self.anchorids[0],1], self.pos[self.npids[0],1],
                                            self.pos[self.anchorids[0],2], self.pos[self.npids[0],2])
                    self.reference = np.array(np.array([(self.R2)**0.5, 0.0, 0.0]) + self.box_length_half)
                    print (self.R2)**0.5
                self.shiftpos[self.npids[0],:] = self.box_length_half
                self.mols[self.npids[0]] = 0
                for igraft in range(self.nanchors):
                    self.mols[igraft] = igraft
                    shift, orientation = self.getGraft(imol, igraft)

                    self.shiftpos[self.anchorids[igraft],:] = self.reference
                    # -orientation = NP center position, 0 is graft
                    for iatom in self.graftids[igraft,:]:
                        self.mols[iatom] = igraft
                        #iatom = self.atomids.index(gatom)
                        #self.shiftpos[iatom,:] = self.pos[iatom,:] - shift[:]
                        # shift to origin-reference
                        self.shiftpos[iatom,:] = self.pos[iatom,:] - self.pos[self.npids[0],:]
                        # un wrap
                        for idim in [0,1,2]:
                            if ( self.shiftpos[iatom,idim] > self.box_length_half[idim] ):
                                self.shiftpos[iatom,idim] -= self.box_length_half[idim]
                            elif ( self.shiftpos[iatom,idim] < -self.box_length_half[idim] ):
                                self.shiftpos[iatom,idim] += self.box_length_half[idim]
                        # rotate
                        self.shiftpos[iatom,:] = orientation.dot( self.shiftpos[iatom,:] )
                        # shift to center
                        self.shiftpos[iatom,:] = self.shiftpos[iatom,:] + self.box_length_half
                        #for idim in [0,1,2]:
                        #    if ( self.shiftpos[iatom,idim] > self.box[idim] ):
                        #        self.shiftpos[iatom,idim] -= self.box_length[idim]
                        #    elif ( self.shiftpos[iatom,idim] < 0.0 ):
                        #        self.shiftpos[iatom,idim] += self.box_length[idim]
                        # orient
                        #self.shiftpos[iatom,:] = orientation.dot( self.shiftpos[iatom,:] )
                        #self.shiftpos[iatom,:] = np.dot(np.dot( orientation, self.shiftpos[iatom,:] ), np.flip(orientation))
                        #self.shiftpos[iatom,:] = np.dot( self.shiftpos[iatom,:], orientation )

            else:
                shift = self.getShift(imol)
                for iatom in self.molids[imol]:
                    self.shiftpos[iatom,:] = self.pos[iatom,:] - shift[:]
                    for idim in [0,1,2]:
                        if ( self.shiftpos[iatom,idim] > self.box[idim, 1] ):
                            self.shiftpos[iatom,idim] -= self.box_length[idim]
                        elif ( self.shiftpos[iatom,idim] < self.box[idim, 0] ):
                            self.shiftpos[iatom,idim] += self.box_length[idim]

    def getShift(self, imol):
        if self.center == 'com':
            return self.getCOM(imol)
        elif self.center == 'end':
            return self.getEnd(imol)
        elif self.center == 'center':
            return self.getCenter(imol)

    def getGraft(self, imol, igraft):
        from scipy.spatial.transform import Rotation as R
        from numpy import linalg as LA
        #shift = self.pos[self.anchorids[igraft], :]
        shift = np.array(self.pos[self.anchorids[igraft],:] - self.reference )

        anchorpoint = self.pos[self.anchorids[igraft], :] 
        vector = np.array(anchorpoint - self.pos[self.npids[imol], :] )
        #vector = np.array(anchorpoint - self.pos[self.npids[imol], :] + self.reference[0])
        nvector = vector/LA.norm(vector)
        ortho = np.cross(nvector, np.array([1,0,0]))
        northo = ortho/LA.norm(ortho)
        theta = np.arccos( np.dot(nvector, np.array([1,0,0])))
        r = np.cos(theta/2.0) + northo * np.sin(theta/2.0)
        #return shift[:], r

        vector = np.array(self.pos[self.anchorids[igraft],:] - self.pos[self.npids[imol], :])
        alpha = np.arccos(-nvector[1]/(1.-nvector[2]**2.)**0.5)
        beta = np.arccos(nvector[2])
        gamma = np.arccos(-nvector[1]/(1.-nvector[2]**2.)**0.5)
        temp = (vector[0]**2.+vector[1]**2.)**0.5
        r = R.from_euler('zx', [np.arcsin(vector[2]/self.reference[0]), np.arccos(temp/self.reference[0])*2. ], degrees=False)
        #return shift[:], r.as_dcm()

        vec1 = np.array(self.pos[self.anchorids[igraft],:] - self.pos[self.npids[imol], :])
        vec2 = np.array([1,0,0])
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        #r = R.from_euler([np.arccos(-vector[1]/(1.-vector[2]**2.)**0.5), 0, np.sin(np.pi/4), np.cos(np.pi/4)])
        
        #R = np.linalg.solve(vector,self.reference).T 
        #R = np.dot( np.flip(vector), self.reference )
        #rotation matrix
        return shift[:], rotation_matrix
    def getCenter(self, imol):
        return self.pos[np.median(self.atomids[imol,:]),:]
    def getEnd(self, imol):
        return self.pos[min(self.atomids[imol,:]),:]
    def getCOM(self, imol):
        com_gnp_pos = np.zeros( (self.NatomsNP[inp], 3), dtype=np.float)
        name = []

        # loop over chains/molecuels
        for iatom in range(self.NatomsNP[inp]):
            for idim in range(3) :
                # shift the positions by the com, so that the com is 0
                com_gnp_pos[idim] = self.pos[iatom,idim] - com[idim]
        return pos

    def comNPPosName(self, imol):
        """
        Calculate the center of mass of a nanoparticle density of 
        all types as a function of NP COM
        returns array of [total, type_1, type_2, ...]
        """
        com_gnp_pos = np.zeros( (self.NatomsNP[inp], 3), dtype=np.float)
        name = []

        # loop over chains/molecuels
        for iatom in range(self.NatomsNP[inp]):
            for idim in range(3) :
                # shift the positions by the com, so that the com is 0
                com_gnp_pos[idim] = self.pos[iatom,idim] - com[idim]
                if ( com_gnp_pos[idim] > self.box[idim, 1] ):
                    com_gnp_pos[idim] -= self.box_length[idim]
                elif ( com_gnp_pos[idim] < self.box[idim, 0] ):
                    com_gnp_pos[idim] += self.box_length[idim]

        return com, com_gnp_pos, name


    def findGrafts(self, frame, imol):
        nps = []
        anchors = []
        grafts = []
        nmol = -1
        for iatom in range(frame['natoms']):
            #if frame['mol_name'][iatom] == self.gnp_mol_name:
            #print frame['atype'][iatom]
            if frame['amol'][iatom] == self.molids[imol]:
                atype = frame['atype'][iatom]
                if str(atype) == self.np_aname:
                    nps.append(iatom)
                elif str(atype) == self.anchor_aname:
                    anchors.append(iatom)
                elif str(atype) == self.graft_aname:
                    grafts.append(iatom)
        self.nanchors = len(anchors)
        # assume the atomids are in order and the only thing between anchors in list are grafts, ie ordered
        sort_anchors = np.sort(anchors)
        graftids = np.zeros( (self.nanchors, len(grafts)), dtype=np.int)
        maxNgraft = 0 
        for i in range(self.nanchors):
            graftfirst = sort_anchors[i] + 1
            if i < self.nanchors-1:
                graftend = sort_anchors[i+1] - 1
            else:
                graftend = max(grafts)
            ianchor = anchors.index(sort_anchors[i])
            igraft = 0
            for j in range(graftfirst,graftend+1):
                graftids[ianchor,igraft] = j
                igraft += 1
            maxNgraft = max(maxNgraft, igraft)
        tgraftids = np.array(graftids[:,:maxNgraft])

        return tgraftids, anchors, nps

    def findCenterIds(self, frame):
        molids = []
        atomids = []
        nmol = -1
        for iatom in range(frame['natoms']):
            #if frame['mol_name'][iatom] == self.gnp_mol_name:
            #print frame['atype'][iatom]
            if frame['atype'][iatom] == self.center_type:
                amol = frame['amol'][iatom]
                if amol not in molids:
                    atomids.append( [] )
                    molids.append( frame['amol'][iatom] )
                    nmol += 1
                atomids[molids.index(amol)-1].append( iatom )

        return atomids, molids

    def findSolutionids(self, frame):
        molids = []
        atomids = []
        nmol = -1
        for iatom in range(frame['natoms']):
            if ((str(frame['atype'][iatom]) in self.solution_anames) or 
               (frame['aname'][iatom] in self.solution_anames)):
                amol = frame['amol'][iatom]
                if amol not in molids:
                    atomids.append( [] )
                    molids.append( amol )
                    nmol += 1
                
                atomids[molids.index(amol)-1].append( iatom )

        return atomids, molids

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Analyze configurations from simulations '
                                                 'of surfactant, chain molecules '
                                                 'which are likely aggregated')
    parser.add_argument('-t', "--traj_file", type=str, nargs='+',
                   help='Configuration/Trajectory/movie file: *xyz '
                        'Right now only 1 file capability, ' 
                        'pdb and towhee_movie are under development')
    parser.add_argument("--traj_type", type=str, choices=['xyz', 'gro', 'pdb','dump'],
                   help='trajectory filetype')
    parser.add_argument("--start_frame", type=int, default=1,
                   help='Set the starting frame number.')
    parser.add_argument("--end_frame", type=int, default=1,
                   help='Set the ending frame number.')
    parser.add_argument('-b', "--box_length", type=float, nargs='+',
                   help='Simulation box length (orthorombic is the only option)'
                        '0-box_lenth can give only 3')
    parser.add_argument("--gnp_name", type=str, default='GP',
                   help='Molecule name for grafted nanoparticle')
    parser.add_argument("--np_name", type=str, default='O',
                   help='charachter name for the NP name')
    parser.add_argument("--anchor_name", type=str, default='X',
                   help='charachter name for the graft anchor name')
    parser.add_argument("--graft_name", type=str, default='H',
                   help='charachter name for the graft name')
    parser.add_argument("--matrix_name", type=str, default='He',
                   help='charachter name for the matrix name')
    parser.add_argument("--center_type", type=int, default=3,
                   help='molecule type you want to center')
    parser.add_argument("-o", "--output", metavar='outputFile', type=str, choices=['chk', 'frag', 'xyz', 'gro', 'g96', 'xml', 'dump'],
                   help='Output a file in a specific format (frag requires mcf to zero at COM): -o xyz')
    parser.add_argument("--center", type=str, choices=['com', 'end', 'center', 'graft'],
                   help='which molecules do you want to set to the center')

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

    #if not analyze_error:
    #    analyze.write()

if __name__ == '__main__':
    sys.exit(main())

