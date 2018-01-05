#!/usr/bin/python
# coding = utf-8
#This module is part of an analysis package

Authorinfo = """
             ------------------Name: Yuan-Chao Hu--------------
             --------------Email: ychu0213@gmail.com-----------
             ----------Web: https://yuanchaohu.github.io/------
             """

Docstr = """
             Reading partcles' positions from snapshots of molecular simulations

             To run the code, dump file format (id type x y z ...) are needed
             It supports three types of coordinates now x, xs, xu
             To obtain the dump information, use 
             ********************* from dump import readdump*********************
             ******************** d = readdump(filename, ndim) ******************
             ************************ d.read_onefile() **************************

             d.TimeStep:        a list of timestep from each snapshot
             d.ParticleNumber:  a list of particle number from each snapshot
             d.ParticleType:    a list of particle type in array in each snapshot
             d.Positions:       a list of particle coordinates in array in each snapshot
             d.Snapshotnumber:  snapshot number 
             d.Boxlength:       a list of box length in array in each snapshot
             d.Boxbounds:       a list of box boundaries in array in each snapshot
             d.hmatrix:         a list of h-matrix of the cells in each snapshot

             The information is stored in list whose elements are mainly numpy arraies

             This module is powerful for different dimensions by giving 'ndim' for orthogonal boxes
             For a triclinic box, convert the bounding box back into the trilinic box parameters:
             xlo = xlo_bound - MIN(0.0,xy,xz,xy+xz)
             xhi = xhi_bound - MAX(0.0,xy,xz,xy+xz)
             ylo = ylo_bound - MIN(0.0,yz)
             yhi = yhi_bound - MAX(0.0,yz)
             zlo = zlo_bound
             zhi = zhi_bound
             See 'http://lammps.sandia.gov/doc/Section_howto.html#howto-12'
         """

import numpy as np 


class readdump:
    """Read snapshots from simulations"""

    def __init__(self, filename, ndim, *arg):
        self.filename       = filename #input snapshots
        self.ndim           = ndim     #dimension 
        self.TimeStep       = []       #simulation timestep @ each snapshot
        self.ParticleNumber = []       #particle's Number @ each snapshot
        self.ParticleType   = []       #particle's type @ each snapshot
        self.Positions      = []       #a list containing all snapshots, each element is a snapshot
        self.SnapshotNumber = 0        #snapshot number
        self.Boxlength      = []       #box length @ each snapshot
        self.Boxbounds      = []       #box boundaries @ each snapshot
        self.Realbounds     = []       #real box bounds of a triclinic box
        self.hmatrix        = []       #h-matrix of the cells

    def read_onefile(self):
        """ Read all snapshots from one dump file """

        f = open(self.filename, 'r')

        positions = self.read_snapshot(f)
        while positions.any() :
            self.Positions.append(positions)
            positions = self.read_snapshot(f)

        f.close()
        self.SnapshotNumber = len(self.TimeStep)
        print ('--------Reading Over---------')

    def read_multiple():
        """ Read all snapshots from individual dump files """
        pass

    def read_snapshot(self, f):
        """ Read a snapshot at one time """

        try:
            item = f.readline()
            timestep = int(f.readline().split()[0])
            self.TimeStep.append(timestep)
            item = f.readline()
            ParticleNumber = int(f.readline())
            self.ParticleNumber.append(ParticleNumber)

            item = f.readline().split()
            #-------Read Orthogonal Boxes---------
            if not 'xy' in item:
                boxbounds = np.zeros((self.ndim, 2))   #box boundaries of (x y z)
                boxlength = np.zeros(self.ndim)   #box length along (x y z)
                for i in range(self.ndim):
                    item = f.readline().split()
                    boxbounds[i, :] = item[:2]

                boxlength = boxbounds[:, 1] - boxbounds[:, 0]           
                if self.ndim < 3:
                    for i in range(3 - self.ndim):
                        f.readline()
                self.Boxbounds.append(boxbounds)
                self.Boxlength.append(boxlength)
                hmatrix = np.diag(boxlength)
                self.hmatrix.append(hmatrix)

                item = f.readline().split()
                names = item[2:]
                positions = np.zeros((ParticleNumber, self.ndim))
                ParticleType = np.zeros(ParticleNumber, dtype=np.int)

                if 'xu' in names: 
                    for i in range(ParticleNumber):
                        item = f.readline().split()
                        ParticleType[int(item[0]) - 1] = int(item[1])
                        positions[int(item[0]) - 1] = [float(j) for j in item[2: self.ndim + 2]]

                elif 'x' in names: 
                    for i in range(ParticleNumber):
                        item = f.readline().split()
                        ParticleType[int(item[0]) - 1] = int(item[1])
                        positions[int(item[0]) - 1] = [float(j) for j in item[2: self.ndim + 2]]

                    positions = np.where(positions < boxbounds[:, 0], positions + boxlength, positions)
                    positions = np.where(positions > boxbounds[:, 1], positions - boxlength, positions)

                elif 'xs' in names: 
                    for i in range(ParticleNumber):
                        item = f.readline().split()
                        ParticleType[int(item[0]) - 1] = int(item[1])
                        positions[int(item[0]) - 1] = [float(j) for j in item[2: self.ndim + 2]] * boxlength

                    positions = np.where(positions < boxbounds[:, 0], positions + boxlength, positions)
                    positions = np.where(positions > boxbounds[:, 1], positions - boxlength, positions)

                self.ParticleType.append(ParticleType)
                return positions

            #-------Read Triclinic Boxes---------
            else:
                boxbounds  = np.zeros((self.ndim, 3))   #box boundaries of (x y z) with tilt factors
                boxlength  = np.zeros(self.ndim)   #box length along (x y z)
                for i in range(self.ndim):
                    item = f.readline().split()
                    boxbounds[i, :] = item[:3]    #with tilt factors
                if self.ndim < 3:
                    for i in range(3 - self.ndim):
                        item = f.readline().split()
                        boxbounds = np.vstack((boxbounds, np.array(item[:3], dtype = np.float)))

                xlo_bound, xhi_bound, xy = boxbounds[0, :]
                ylo_bound, yhi_bound, xz = boxbounds[1, :]
                zlo_bound, zhi_bound, yz = boxbounds[2, :]
                xlo = xlo_bound - min((0.0, xy, xz, xy + xz))
                xhi = xhi_bound - max((0.0, xy, xz, xy + xz))
                ylo = ylo_bound - min((0.0, yz))
                yhi = yhi_bound - max((0.0, yz))
                zlo = zlo_bound
                zhi = zhi_bound
                h0  = xhi - xlo
                h1  = yhi - ylo
                h2  = zhi - zlo 
                h3  = yz 
                h4  = xz 
                h5  = xy 

                realbounds = np.array([xlo, xhi, ylo, yhi, zlo, zhi]).reshape((3, 2))
                self.Realbounds.append(realbounds[:self.ndim])
                reallength = (realbounds[:, 1] - realbounds[:, 0])[:self.ndim]
                self.Boxlength.append(reallength)
                boxbounds  = boxbounds[:self.ndim, :2]
                self.Boxbounds.append(boxbounds)
                hmatrix = np.zeros((3, 3))
                hmatrix[0] = [h0, 0 , 0]
                hmatrix[1] = [h5, h1, 0]
                hmatrix[2] = [h4, h3, h2]
                hmatrix    = hmatrix[:self.ndim]
                self.hmatrix.append(hmatrix)

                item = f.readline().split()
                names = item[2:]
                positions = np.zeros((ParticleNumber, self.ndim))
                ParticleType = np.zeros(ParticleNumber, dtype=np.int)
                if 'x' in names:
                    for i in range(ParticleNumber):
                        item = f.readline().split()
                        ParticleType[int(item[0]) - 1] = int(item[1])
                        positions[int(item[0]) - 1] = [float(j) for j in item[2: self.ndim + 2]]

                elif 'xs' in names: 
                    for i in range(ParticleNumber):
                        item = f.readline().split()
                        pid  = int(item[0]) - 1
                        ParticleType[pid] = int(item[1])
                        if self.ndim == 3:
                            positions[pid, 0] = xlo_bound + float(item[2])*h0 + float(item[3])*h5 + float(item[4])*h4
                            positions[pid, 1] = ylo_bound + float(item[3])*h1 + float(item[4])*h3
                            positions[pid, 2] = zlo_bound + float(item[4])*h2
                        elif self.ndim == 2:
                            positions[pid, 0] = xlo_bound + float(item[2])*h0 + float(item[3])*h5 
                            positions[pid, 1] = ylo_bound + float(item[3])*h1

                self.ParticleType.append(ParticleType)
                return positions

        except:
            positions= np.zeros((self.ParticleNumber[0], self.ndim))
            return positions
