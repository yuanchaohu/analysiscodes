# coding = utf8
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
             ********************** d = readdump(filename) **********************
             ************************ d.read_onefile() **************************

             d.TimeStep:        a list of timestep from each snapshot
             d.ParticleNumber:  a list of particle number from each snapshot
             d.ParticleType:    a list of particle type in array in each snapshot
             d.Positions:       a list of particle coordinates in array in each snapshot
             d.Snapshotnumber:  snapshot number 
             d.Boxlength:       a list of box length in array in each snapshot
             d.Boxbounds:       a list of box boundaries in array in each snapshot

             The information is stored in list whose element is mainly numpy array

             This code is usable for both 3D and 2D cases, but in the latter case
             (x, y, z) coordinates are all needed
         """

import numpy as np 

class readdump:
    """Read snapshots of a simulation"""

    def __init__(self, filename, *arg):
        self.filename = filename
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

            boxbounds = np.zeros((3, 2), dtype=np.float)   #box boundaries of (x y z)
            boxlength = np.zeros(3, dtype=np.float)   #box length along (x y z)

            item = f.readline().split()
            boxbounds[0, :] = item[:2]
            boxlength[0] = float(item[1]) - float(item[0])

            item = f.readline().split()
            boxbounds[1, :] = item[:2]
            boxlength[1] = float(item[1]) - float(item[0])

            item = f.readline().split()
            boxbounds[2, :] = item[:2]
            boxlength[2] = float(item[1]) - float(item[0])

            self.Boxbounds.append(boxbounds)
            self.Boxlength.append(boxlength)

            item = f.readline().split()
            names = item[2:]
            positions = np.zeros((ParticleNumber, 3), dtype=np.float)
            ParticleType = np.zeros(ParticleNumber, dtype=np.int)

            if 'x' in names: 
                #print ('----------Reading unscaled coordinates x---------')
                for i in range(ParticleNumber):
                    item = f.readline().split()
                    ParticleType[int(item[0]) - 1] = int(item[1])
                    positions[int(item[0]) - 1, 0] = float(item[2])
                    positions[int(item[0]) - 1, 1] = float(item[3])
                    positions[int(item[0]) - 1, 2] = float(item[4])

            elif 'xu' in names: 
                #print ('----------Reading unscaled coordinates xu---------')
                for i in range(ParticleNumber):
                    item = f.readline().split()
                    ParticleType[int(item[0]) - 1] = int(item[1])
                    positions[int(item[0]) - 1, 0] = float(item[2])
                    positions[int(item[0]) - 1, 1] = float(item[3])
                    positions[int(item[0]) - 1, 2] = float(item[4])

            elif 'xs' in names: 
                #print ('----------Reading scaled coordinates xs---------')
                for i in range(ParticleNumber):
                    item = f.readline().split()
                    ParticleType[int(item[0]) - 1] = int(item[1])
                    positions[int(item[0]) - 1, 0] = float(item[2]) * float(boxlength[0])
                    positions[int(item[0]) - 1, 1] = float(item[3]) * float(boxlength[1])
                    positions[int(item[0]) - 1, 2] = float(item[4]) * float(boxlength[2])

            self.ParticleType.append(ParticleType)
            return  positions

        except:
            positions= np.zeros((self.ParticleNumber[0], 3), dtype=np.float)
            return positions