#!/usr/bin/python
# coding = utf-8
#This module is part of an analysis package

Authorinfo = """
             ----------Name: Qi Liu [lq4866@gmail.com]--------------------
             ----------Name: Yuan-Chao Hu [ychu0213@gmail.com]------------
             ----------Web: https://yuanchaohu.github.io/-----------------
             """

Docstr = """
         This module is responsible for calculating pair entropy S2
         The code is written refering to 'Allen book P183' covering unary to Quinary systems

         The code accounts for both orthogonal and triclinic cells
         """

import os
import numpy as np 
from   dump  import readdump

def Nidealfac(ndim):
    """ Choose factor of Nideal in g(r) calculation """

    if ndim == 3:
        return 4.0 / 3
    elif ndim == 2:
        return 1.0

def Areafac(ndim):
    """ Choose factor of area in S2 calculation """

    if ndim == 3:
        return 4.0  #4 * PI * R2
    elif ndim == 2:
        return 2.0  #2 * PI * R

class S2:
    """ Compute pair entropy S2 """

    def __init__(self, inputfile, ndim, *arg):
        self.inputfile = inputfile
        self.ndim = ndim 
        d = readdump(self.inputfile, self.ndim)
        d.read_onefile()

        self.TimeStep     = d.TimeStep[1] - d.TimeStep[0]
        if self.TimeStep != d.TimeStep[-1] - d.TimeStep[-2]:
            raise ValueError('*********** dump interval changes **************') 
        self.ParticleNumber     = d.ParticleNumber[0] 
        if d.ParticleNumber[0] != d.ParticleNumber[-1]:
            raise ValueError('************* Paticle Number Changes **************')
        self.ParticleType   = d.ParticleType
        self.Positions      = d.Positions
        self.SnapshotNumber = d.SnapshotNumber
        self.Boxlength      = d.Boxlength[0]
        if not (d.Boxlength[0] == d.Boxlength[-1]).all():
            raise ValueError('*********Box Length Changed from Dump***********')
        self.hmatrix    = d.hmatrix
        self.Volume     = np.prod(self.Boxlength)
        self.rhototal   = self.ParticleNumber / self.Volume #number density of the total system
        self.typecounts = np.unique(self.ParticleType[0], return_counts = True) 
        self.Type       = self.typecounts[0]
        self.TypeNumber = self.typecounts[1]
        self.rhotype    = self.TypeNumber / self.Volume
        print ('Particle Type:', self.Type)
        print ('Particle TypeNumber:', self.TypeNumber)
        if np.sum(self.TypeNumber) != self.ParticleNumber:
            raise ValueError('****** Sum of Indivdual Types is Not the Total Amount*******')

    def timeave(self, outputfile, ppp, avetime, rdelta = 0.01, dt = 0.002, results_path='../../analysis/S2/'):
        """ Calculate Time Averaged S2
            
            outputfile is file name storing outputs
            ppp is periodic boundary conditions, set 1 for yes and 0 for no, should be a list
            avetime is the time scale used to average S2, in time unit; set 0 to get results of individual snapshot
            rdelta is the bin size of gr, default is 0.01
            dt is timestep of a simulation, default is 0.002
            results_path is the path of outputfile, default is '../../analysis/S2/'
        """
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        if len(self.Type) == 1:
            S2results      = self.Unary(ppp, rdelta)
        if len(self.Type) == 2: 
            S2results      = self.Binary(ppp, rdelta)
        if len(self.Type) == 3: 
            S2results =    self.Ternary(ppp, rdelta)
        if len(self.Type) == 4: 
            S2results      = self.Quarternary(ppp, rdelta)
        if len(self.Type) == 5: 
            S2results      = self.Quinary(ppp, rdelta)

        if avetime:
            avetime    = int(avetime / dt / self.TimeStep)
            S2average  = np.zeros((self.ParticleNumber, self.SnapshotNumber - avetime)) 
            for n in range(self.SnapshotNumber - avetime):
                S2average[:, n] = S2results[:, n: n + avetime].mean(axis = 1)
            results = np.column_stack((np.arange(self.ParticleNumber) + 1, S2average))
        else:
            results = np.column_stack((np.arange(self.ParticleNumber) + 1, S2results))

        names      = 'id   S2_of_each_snapshot'
        fileformat = '%d ' + '%.6f ' * (self.SnapshotNumber - avetime)
        np.savetxt(results_path + outputfile, results, fmt = fileformat, header = names, comments = '')
        print ('---------- Get S2 results over ---------')
        return results

    def spatialcorr(self, outputcorr, outputs2, ppp, avetime, rdelta = 0.01, dt = 0.002, results_path='../../analysis/S2/'):
        """ Calculate Spatial Correlation of Time Averaged S2
            
            Excuting this function will excute the function timeave() first, so the averaged S2 will be output
            outputcorr is file name storing outputs of spatial correlation of time averaged S2
            outputs2 is file name storing outputs of time averaged S2
            ppp is periodic boundary conditions, set 1 for yes and 0 for no, should be a list
            avetime is the time scale used to average S2, in time unit; set 0 to get results of individual snapshot
            rdelta is the bin size of gr, default is 0.01
            dt is timestep of a simulation, default is 0.002
            results_path is the path of outputfile, default is '../../analysis/S2/'
        """
        print ('--------- Calculating Spatial Correlation of S2 ------')
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        timeaves2      = self.timeave(outputs2, ppp, avetime, rdelta, dt, results_path)[:, 1:].T 
        MAXBIN         = int(self.Boxlength.min() / 2.0 / rdelta)
        grresults      = np.zeros((MAXBIN, 3))
        SnapshotNumber = len(timeaves2[:, 0])
        Positions      = self.Positions[:SnapshotNumber]
        for n in range(SnapshotNumber):
            hmatrixinv = np.linalg.inv(self.hmatrix[n])
            for i in range(self.ParticleNumber - 1):
                RIJ      = Positions[n][i+1:] - Positions[n][i]
                #periodic = np.where(np.abs(RIJ / self.Boxlength[np.newaxis, :]) > 0.50, np.sign(RIJ), 0).astype(np.int)
                #RIJ -= self.Boxlength * periodic * ppp    #remove PBC
                matrixij = np.dot(RIJ, hmatrixinv)
                RIJ      = np.dot(matrixij - np.rint(matrixij) * ppp, self.hmatrix[n]) #remove PBC
                distance = np.sqrt(np.square(RIJ).sum(axis = 1))
                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 0]    += Countvalue
                S2IJ                = timeaves2[n, i + 1:] * timeaves2[n, i]
                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta), weights = S2IJ)
                grresults[:, 1]    += Countvalue

        binleft          = BinEdge[:-1]   #real value of each bin edge, not index
        binright         = BinEdge[1:]    #len(Countvalue) = len(BinEdge) - 1
        Nideal           = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim) 
        grresults       *= 2 / self.ParticleNumber / SnapshotNumber / (Nideal[:, np.newaxis] * self.rhototal)
        grresults[:, 2]  = np.where(grresults[:, 0] != 0, grresults[:, 1] / grresults[:, 0], np.nan)
        binright         = binright - 0.5 * rdelta #middle of each bin
        results          = np.column_stack((binright, grresults))
        names            = 'r   g(r)   gs2(r)   gs2/gr'
        np.savetxt(results_path + outputcorr, results, fmt='%.6f', header = names, comments = '')
        print ('---------- Get gs2(r) results over ---------')
        return results

    def Unary(self, ppp, rdelta):
        print ('--------- This is a Unary System ------')

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        S2results   = np.zeros((self.ParticleNumber, self.SnapshotNumber))
        for n in range(self.SnapshotNumber):
            hmatrixinv = np.linalg.inv(self.hmatrix[n])
            for i in range(self.ParticleNumber):
                RIJ          = np.delete(self.Positions[n], i, axis = 0) - self.Positions[n][i]
                matrixij     = np.dot(RIJ, hmatrixinv)
                RIJ          = np.dot(matrixij - np.rint(matrixij) * ppp, self.hmatrix[n]) #remove periodic boundary conditions
                distance     = np.sqrt(np.square(RIJ).sum(axis = 1))
                
                particlegr          = np.zeros(MAXBIN) #11
                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta)) #1-1; 
                binleft             = BinEdge[:-1]          #real value of each bin edge, not index
                binright            = BinEdge[1:]           #len(Countvalue) = len(BinEdge) - 1
                Nideal              = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
                particlegr          = Countvalue / Nideal / self.rhototal
                integralgr          = (particlegr * np.log(particlegr + 1e-12) - (particlegr - 1)) * self.rhototal
                binright           -= 0.5 * rdelta                                  #middle of each bin
                S2results[i, n]     =-0.5 * np.sum(Areafac(self.ndim) * np.pi * binright**(self.ndim - 1) * integralgr * rdelta)

        return S2results

    def Binary(self, ppp, rdelta):
        print ('--------- This is a Binary System ------')

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        S2results   = np.zeros((self.ParticleNumber, self.SnapshotNumber))
        for n in range(self.SnapshotNumber):
            hmatrixinv = np.linalg.inv(self.hmatrix[n])
            for i in range(self.ParticleNumber):
                RIJ          = np.delete(self.Positions[n], i, axis = 0) - self.Positions[n][i]
                matrixij     = np.dot(RIJ, hmatrixinv)
                RIJ          = np.dot(matrixij - np.rint(matrixij) * ppp, self.hmatrix[n]) #remove periodic boundary conditions
                distance     = np.sqrt(np.square(RIJ).sum(axis = 1))
                particletype = np.delete(self.ParticleType[n], i)
                TIJ          = np.c_[particletype, np.zeros_like(particletype) + self.ParticleType[n][i]]
                Countsum     = TIJ.sum(axis = 1)
                
                particlegr          = np.zeros((MAXBIN, 3)) #11 12/21 22
                Countvalue, BinEdge = np.histogram(distance[Countsum == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta)) #1-1; 
                binleft             = BinEdge[:-1]          #real value of each bin edge, not index
                binright            = BinEdge[1:]           #len(Countvalue) = len(BinEdge) - 1
                Nideal              = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
                particlegr[:, 0]    = Countvalue / Nideal / self.rhotype[0]
                Countvalue, BinEdge = np.histogram(distance[Countsum == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta)) #1-2; 2-1
                rho12               = self.Type[self.Type != self.ParticleType[n][i]] - 1
                particlegr[:, 1]    = Countvalue / Nideal / self.rhotype[rho12]
                Countvalue, BinEdge = np.histogram(distance[Countsum == 4], bins = MAXBIN, range = (0, MAXBIN * rdelta)) #2-2
                particlegr[:, 2]    = Countvalue / Nideal / self.rhotype[1]
                integralgr          = (particlegr * np.log(particlegr + 1e-12) - (particlegr - 1)) * [self.rhotype[0], self.rhotype[rho12], self.rhotype[1]]
                integralgr          = integralgr[:, np.any(particlegr, axis = 0)]   #remove zero columns in gr 
                binright           -= 0.5 * rdelta                                  #middle of each bin
                S2results[i, n]     =-0.5 * np.sum(Areafac(self.ndim) * np.pi * binright**(self.ndim - 1) * integralgr.sum(axis = 1) * rdelta)

        return S2results

    def Ternary(self, ppp, rdelta):
        print ('--------- This is a Ternary System ------')

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        S2results   = np.zeros((self.ParticleNumber, self.SnapshotNumber))
        for n in range(self.SnapshotNumber):
            hmatrixinv = np.linalg.inv(self.hmatrix[n])
            for i in range(self.ParticleNumber):
                RIJ          = np.delete(self.Positions[n], i, axis = 0) - self.Positions[n][i]
                matrixij     = np.dot(RIJ, hmatrixinv)
                RIJ          = np.dot(matrixij - np.rint(matrixij) * ppp, self.hmatrix[n]) #remove periodic boundary conditions
                distance     = np.sqrt(np.square(RIJ).sum(axis = 1))
                particletype = np.delete(self.ParticleType[n], i)
                TIJ          = np.c_[particletype, np.zeros_like(particletype) + self.ParticleType[n][i]]
                Countsum     = TIJ.sum(axis = 1)
                Countsub     = np.abs(TIJ[:, 0] - TIJ[:, 1])

                particlegr          = np.zeros((MAXBIN, 6)) #11 12/21 13/31 22 23/32 33
                usedrho             = np.zeros(6)
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                binleft             = BinEdge[:-1]          #real value of each bin edge, not index
                binright            = BinEdge[1:]           #len(Countvalue) = len(BinEdge) - 1
                Nideal              = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
                particlegr[:, 0]    = Countvalue / Nideal / self.rhotype[0] #11
                usedrho[0]          = self.rhotype[0]
                Countvalue, BinEdge = np.histogram(distance[(Countsum == 4) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 1]    = Countvalue / Nideal / self.rhotype[1] #22
                usedrho[1]          = self.rhotype[1]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 6], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 2]    = Countvalue / Nideal / self.rhotype[2] #33
                usedrho[2]          = self.rhotype[2]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho12               = self.Type[self.Type != self.ParticleType[n][i]][0] - 1
                particlegr[:, 3]    = Countvalue / Nideal / self.rhotype[rho12] #12/21
                usedrho[3]          = self.rhotype[rho12]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 5], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho23               = self.Type[self.Type != self.ParticleType[n][i]][1] - 1
                particlegr[:, 4]    = Countvalue / Nideal / self.rhotype[rho23] #23/32
                usedrho[4]          = self.rhotype[rho23]
                Countvalue, BinEdge = np.histogram(distance[(Countsum == 4) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]]
                rho13               = medium[medium != 2] - 1
                particlegr[:, 5]    = Countvalue / Nideal / self.rhotype[rho13] #13/31
                usedrho[5]          = self.rhotype[rho13]

                integralgr          = (particlegr * np.log(particlegr + 1e-12) - (particlegr - 1)) * usedrho
                integralgr          = integralgr[:, np.any(particlegr, axis = 0)]   #remove zero columns in gr 
                binright           -= 0.5 * rdelta                                  #middle of each bin
                S2results[i, n]     =-0.5 * np.sum(Areafac(self.ndim) * np.pi * binright**(self.ndim - 1) * integralgr.sum(axis = 1) * rdelta)

        return S2results

    def Quarternary(self, ppp, rdelta):
        print ('--------- This is a Quarternary System ------')

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        S2results   = np.zeros((self.ParticleNumber, self.SnapshotNumber))
        for n in range(self.SnapshotNumber):
            hmatrixinv = np.linalg.inv(self.hmatrix[n])
            for i in range(self.ParticleNumber):
                RIJ          = np.delete(self.Positions[n], i, axis = 0) - self.Positions[n][i]
                matrixij     = np.dot(RIJ, hmatrixinv)
                RIJ          = np.dot(matrixij - np.rint(matrixij) * ppp, self.hmatrix[n]) #remove periodic boundary conditions
                distance     = np.sqrt(np.square(RIJ).sum(axis = 1))
                particletype = np.delete(self.ParticleType[n], i)
                TIJ          = np.c_[particletype, np.zeros_like(particletype) + self.ParticleType[n][i]]
                Countsum     = TIJ.sum(axis = 1)
                Countsub     = np.abs(TIJ[:, 0] - TIJ[:, 1])

                particlegr          = np.zeros((MAXBIN, 10)) 
                usedrho             = np.zeros(10)
                Countvalue, BinEdge = np.histogram(distance[Countsum == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                binleft             = BinEdge[:-1]          #real value of each bin edge, not index
                binright            = BinEdge[1:]           #len(Countvalue) = len(BinEdge) - 1
                Nideal              = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
                particlegr[:, 0]    = Countvalue / Nideal / self.rhotype[0] #11
                usedrho[0]          = self.rhotype[0]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 1]    = Countvalue / Nideal / self.rhotype[1] #22
                usedrho[1]          = self.rhotype[1]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 2]    = Countvalue / Nideal / self.rhotype[2] #33
                usedrho[2]          = self.rhotype[2]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 8], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 3]    = Countvalue / Nideal / self.rhotype[3] #44
                usedrho[3]          = self.rhotype[3]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho12               = self.Type[self.Type != self.ParticleType[n][i]][0] - 1
                particlegr[:, 4]    = Countvalue / Nideal / self.rhotype[rho12] #12
                usedrho[4]          = self.rhotype[rho12]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]][:2]
                rho13               = medium[medium != 2] - 1
                particlegr[:, 5]    = Countvalue / Nideal / self.rhotype[rho13] #13
                usedrho[5]          = self.rhotype[rho13]
                Countvalue, BinEdge = np.histogram(distance[Countsub  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]]
                rho14               = medium[(medium != 2) & (medium != 3)] - 1
                particlegr[:, 6]    = Countvalue / Nideal / self.rhotype[rho14] #14
                usedrho[6]          = self.rhotype[rho14]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 5) & (Countsub == 1)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho23               = self.Type[self.Type != self.ParticleType[n][i]][1] - 1
                particlegr[:, 7]    = Countvalue / Nideal / self.rhotype[rho23] #23
                usedrho[7]          = self.rhotype[rho23]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]][1:]
                rho24               = medium[medium != 3] - 1
                particlegr[:, 8]    = Countvalue / Nideal / self.rhotype[rho24] #24
                usedrho[8]          = self.rhotype[rho24]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 7], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho34               = self.Type[self.Type != self.ParticleType[n][i]][-1] - 1
                particlegr[:, 9]    = Countvalue / Nideal / self.rhotype[rho34] #34
                usedrho[9]          = self.rhotype[rho34]

                integralgr          = (particlegr * np.log(particlegr + 1e-12) - (particlegr - 1)) * usedrho
                integralgr          = integralgr[:, np.any(particlegr, axis = 0)]   #remove zero columns in gr 
                binright           -= 0.5 * rdelta                                  #middle of each bin
                S2results[i, n]     =-0.5 * np.sum(Areafac(self.ndim) * np.pi * binright**(self.ndim - 1) * integralgr.sum(axis = 1) * rdelta)

        return S2results

    def Quinary(self, ppp, rdelta):
        print ('--------- This is a Quinary System ------')

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        S2results   = np.zeros((self.ParticleNumber, self.SnapshotNumber))
        for n in range(self.SnapshotNumber):
            hmatrixinv = np.linalg.inv(self.hmatrix[n])
            for i in range(self.ParticleNumber):
                RIJ          = np.delete(self.Positions[n], i, axis = 0) - self.Positions[n][i]
                matrixij     = np.dot(RIJ, hmatrixinv)
                RIJ          = np.dot(matrixij - np.rint(matrixij) * ppp, self.hmatrix[n]) #remove periodic boundary conditions
                distance     = np.sqrt(np.square(RIJ).sum(axis = 1))
                particletype = np.delete(self.ParticleType[n], i)
                TIJ          = np.c_[particletype, np.zeros_like(particletype) + self.ParticleType[n][i]]
                Countsum     = TIJ.sum(axis = 1)
                Countsub     = np.abs(TIJ[:, 0] - TIJ[:, 1])

                particlegr          = np.zeros((MAXBIN, 15)) 
                usedrho             = np.zeros(15)
                Countvalue, BinEdge = np.histogram(distance[Countsum == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                binleft             = BinEdge[:-1]          #real value of each bin edge, not index
                binright            = BinEdge[1:]           #len(Countvalue) = len(BinEdge) - 1
                Nideal              = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
                particlegr[:, 0]    = Countvalue / Nideal / self.rhotype[0] #11
                usedrho[0]          = self.rhotype[0]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 1]    = Countvalue / Nideal / self.rhotype[1] #22
                usedrho[1]          = self.rhotype[1]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 2]    = Countvalue / Nideal / self.rhotype[2] #33
                usedrho[2]          = self.rhotype[2]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 8) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 3]    = Countvalue / Nideal / self.rhotype[3] #44
                usedrho[3]          = self.rhotype[3]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 10], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                particlegr[:, 4]    = Countvalue / Nideal / self.rhotype[4] #55
                usedrho[4]          = self.rhotype[4]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho12               = self.Type[self.Type != self.ParticleType[n][i]][0] - 1
                particlegr[:, 5]    = Countvalue / Nideal / self.rhotype[rho12] #12
                usedrho[5]          = self.rhotype[rho12]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]][:2]
                rho13               = medium[medium != 2] - 1
                particlegr[:, 6]    = Countvalue / Nideal / self.rhotype[rho13] #13
                usedrho[6]          = self.rhotype[rho13]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 5) & (Countsub == 3)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]][:3]
                rho14               = medium[(medium != 2) & (medium != 3)] - 1
                particlegr[:, 7]    = Countvalue / Nideal / self.rhotype[rho14] #14
                usedrho[7]          = self.rhotype[rho14]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 4)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]]
                rho15               = medium[(medium != 2) & (medium != 3) & (medium != 4)] - 1
                particlegr[:, 8]    = Countvalue / Nideal / self.rhotype[rho15] #15
                usedrho[8]          = self.rhotype[rho15]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 5) & (Countsub == 1)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho23               = self.Type[self.Type != self.ParticleType[n][i]][1] - 1
                particlegr[:,9]     = Countvalue / Nideal / self.rhotype[rho23] #23
                usedrho[9]          = self.rhotype[rho23]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]][1:3]
                rho24               = medium[medium != 3] - 1
                particlegr[:,10]    = Countvalue / Nideal / self.rhotype[rho24] #24
                usedrho[10]         = self.rhotype[rho24]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 7) & (Countsub == 3)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]][1:]
                rho25               = medium[(medium != 3) & (medium != 4)] - 1
                particlegr[:,11]    = Countvalue / Nideal / self.rhotype[rho25] #25
                usedrho[11]         = self.rhotype[rho25]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 7) & (Countsub == 1)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho34               = self.Type[self.Type != self.ParticleType[n][i]][2] - 1
                particlegr[:,12]    = Countvalue / Nideal / self.rhotype[rho34] #34
                usedrho[12]         = self.rhotype[rho34]
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 8) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                medium              = self.Type[self.Type != self.ParticleType[n][i]][2:]
                rho35               = medium[medium != 4] - 1
                particlegr[:,13]    = Countvalue / Nideal / self.rhotype[rho35] #35
                usedrho[13]         = self.rhotype[rho35]
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 9], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                rho45               = self.Type[self.Type != self.ParticleType[n][i]][-1] - 1
                particlegr[:,14]    = Countvalue / Nideal / self.rhotype[rho45] #45
                usedrho[14]         = self.rhotype[rho45]

                integralgr          = (particlegr * np.log(particlegr + 1e-12) - (particlegr - 1)) * usedrho
                integralgr          = integralgr[:, np.any(particlegr, axis = 0)]   #remove zero columns in gr 
                binright           -= 0.5 * rdelta                                  #middle of each bin
                S2results[i, n]     =-0.5 * np.sum(Areafac(self.ndim) * np.pi * binright**(self.ndim - 1) * integralgr.sum(axis = 1) * rdelta)

        return S2results

        