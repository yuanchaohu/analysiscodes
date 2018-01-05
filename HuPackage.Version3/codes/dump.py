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

             The information is stored in list whose elements are mainly numpy arraies

             This module is powerful for different dimensions by giving 'ndim'
         """

import numpy as np 

class readdump:
    """Read snapshots of a simulation"""

    def __init__(self, filename, ndim, *arg):
        self.filename = filename
        self.ndim = ndim            #dimension 
        self.TimeStep = []          #simulation timestep @ each snapshot
        self.ParticleNumber = []    #particle's Number @ each snapshot
        self.ParticleType = []      #particle's type @ each snapshot
        self.Positions = []         #a list containing all snapshots, each element is a snapshot
        self.SnapshotNumber = 0     #snapshot number
        self.Boxlength = []         #box length @ each snapshot
        self.Boxbounds = []         #box boundaries @ each snapshot

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
            item = f.readline()

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

            item = f.readline().split()
            names = item[2:]
            positions = np.zeros((ParticleNumber, self.ndim))
            ParticleType = np.zeros(ParticleNumber, dtype=np.int)

            if ('x' in names) or ('xu' in names): 
                for i in range(ParticleNumber):
                    item = f.readline().split()
                    ParticleType[int(item[0]) - 1] = int(item[1])
                    positions[int(item[0]) - 1] = [float(j) for j in item[2: self.ndim + 2]]

            elif 'xs' in names: 
                for i in range(ParticleNumber):
                    item = f.readline().split()
                    ParticleType[int(item[0]) - 1] = int(item[1])
                    positions[int(item[0]) - 1] = [float(j) for j in item[2: self.ndim + 2]] * boxlength

            self.ParticleType.append(ParticleType)
            return  positions

        except:
            positions= np.zeros((self.ParticleNumber[0], self.ndim))
            return positions
