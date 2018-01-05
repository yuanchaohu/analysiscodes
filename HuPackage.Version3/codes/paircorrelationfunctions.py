#!/usr/bin/python
# coding = utf-8
#This module is part of an analysis package

Authorinfo = """
             ------------------Name: Yuan-Chao Hu--------------
             --------------Email: ychu0213@gmail.com-----------
             ----------Web: https://yuanchaohu.github.io/------
             """

Docstr = """
         This module is responsible for calculating pair correlation functions
         The code is written refering to 'Allen book P183' covering unary to senary systems
         Use the code by 'from paircorrelationfunctions import gr'
         gr(inputfile, ndim).getresults(outputfile, rdelta, ppp, results_path)
         Then the overall and partial g(r) will be written to a file
         """

import os, sys
import numpy as np 
from dump import readdump

def Nidealfac(ndim):
    """ Choose factor of Nideal in g(r) calculation """

    if ndim == 3:
        return 4.0 / 3
    elif ndim == 2:
        return 1.0 

class gr:
    """ Compute (Partial) pair correlation functions g(r) """

    def __init__(self, inputfile, ndim, *arg):
        self.inputfile = inputfile
        self.ndim = ndim 
        d = readdump(self.inputfile, self.ndim)
        d.read_onefile()

        self.ParticleNumber     = d.ParticleNumber[0] 
        if d.ParticleNumber[0] != d.ParticleNumber[-1]:
            raise ValueError('************* Paticle Number Changes **************')
        self.ParticleType   = d.ParticleType
        self.Positions      = d.Positions
        self.SnapshotNumber = d.SnapshotNumber
        self.Boxlength      = d.Boxlength[0]
        if not (d.Boxlength[0] == d.Boxlength[-1]).all():
            raise ValueError('*********Box Length Changed from Dump***********')
        self.Boxbounds  = d.Boxbounds[0]
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

    def getresults(self, outputfile, ppp, rdelta = 0.01, results_path='../../analysis/gr/'):

        if len(self.Type) == 1:
            self.Unary(outputfile, ppp, rdelta, results_path)

        if len(self.Type) == 2: 
            self.Binary(outputfile, ppp, rdelta, results_path)

        if len(self.Type) == 3: 
            self.Ternary(outputfile, ppp, rdelta, results_path)

        if len(self.Type) == 4: 
            self.Quarternary(outputfile, ppp, rdelta, results_path)

        if len(self.Type) == 5: 
            self.Quinary(outputfile, ppp, rdelta, results_path)

        if len(self.Type) == 6: 
            self.Senary(outputfile, ppp, rdelta, results_path)

    def Unary(self, outputfile, ppp, rdelta = 0.01, results_path='../../analysis/gr/'):
        print ('--------- This is a Unary System ------')
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        grresults   = np.zeros(MAXBIN)
        for n in range(self.SnapshotNumber):
            for i in range(self.ParticleNumber - 1):
                RIJ      = self.Positions[n][i+1:] - self.Positions[n][i]
                periodic = np.where(np.abs(RIJ / self.Boxlength) > 0.50, np.sign(RIJ), 0).astype(np.int)
                RIJ -= self.Boxlength * periodic * ppp    #remove PBC
                #RIJ -= np.dot(np.round(RIJ / self.Boxlength[np.newaxis, :]), np.diag(self.Boxlength)) * ppp #another way to remove PBC
                distance = np.sqrt(np.square(RIJ).sum(axis = 1))

                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults += Countvalue

        binleft    = BinEdge[:-1]   #real value of each bin edge, not index
        binright   = BinEdge[1:]   #len(Countvalue) = len(BinEdge) - 1
        Nideal     = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim) 
        grresults  = grresults * 2 / self.ParticleNumber / self.SnapshotNumber / (Nideal * self.rhototal)

        binright = binright - 0.5 * rdelta #middle of each bin
        results  = np.column_stack((binright, grresults))
        names    = 'r  g(r)'
        np.savetxt(results_path + outputfile, results, fmt='%.6f', header = names, comments = '')
        print ('---------- Get g(r) results over ---------')
        return results

    def Binary(self, outputfile, ppp, rdelta = 0.01, results_path='../../analysis/gr/'):
        print ('--------- This is a Binary System ------')
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        grresults   = np.zeros((MAXBIN, 4))
        for n in range(self.SnapshotNumber):
            for i in range(self.ParticleNumber - 1):
                RIJ = self.Positions[n][i+1:] - self.Positions[n][i]
                TIJ = np.c_[self.ParticleType[n][i+1:], np.zeros_like(self.ParticleType[n][i+1:]) + self.ParticleType[n][i]]
                periodic = np.where(np.abs(RIJ / self.Boxlength) > 0.50, np.sign(RIJ), 0).astype(np.int)
                RIJ -= self.Boxlength * periodic * ppp    #remove PBC
                distance = np.sqrt(np.square(RIJ).sum(axis = 1))

                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 0] += Countvalue

                Countsum = TIJ.sum(axis = 1)
                Countvalue, BinEdge = np.histogram(distance[Countsum == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 1] += Countvalue
                Countvalue, BinEdge = np.histogram(distance[Countsum == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 2] += Countvalue
                Countvalue, BinEdge = np.histogram(distance[Countsum == 4], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 3] += Countvalue

        binleft  = BinEdge[:-1]   #real value of each bin edge, not index
        binright = BinEdge[1:]   #len(Countvalue) = len(BinEdge) - 1
        print (self.ndim)
        Nideal   = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
        grresults[:, 0]  = grresults[:, 0] * 2 / self.SnapshotNumber/ self.ParticleNumber / (Nideal * self.rhototal)
        grresults[:, 1]  = grresults[:, 1] * 2 / self.SnapshotNumber/ self.TypeNumber[0] / (Nideal * self.rhotype[0])
        grresults[:, 2]  = grresults[:, 2] * 2 / self.SnapshotNumber/ Nideal * self.Volume / self.TypeNumber[0]/ self.TypeNumber[1]/ 2.0
        grresults[:, 3]  = grresults[:, 3] * 2 / self.SnapshotNumber/ self.TypeNumber[1] / (Nideal * self.rhotype[1])

        binright = binright - 0.5 * rdelta #middle of each bin
        results  = np.column_stack((binright, grresults))
        names    = 'r  g(r)  g11(r)  g12(r)  g22(r)'
        np.savetxt(results_path + outputfile, results, fmt='%.6f', header = names, comments = '')
        print ('---------- Get g(r) results over ---------')
        return results

    def Ternary(self, outputfile, ppp, rdelta = 0.01, results_path='../../analysis/gr/'):
        print ('--------- This is a Ternary System ------')
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        grresults   = np.zeros((MAXBIN, 7))
        for n in range(self.SnapshotNumber):
            for i in range(self.ParticleNumber - 1):
                RIJ = self.Positions[n][i+1:] - self.Positions[n][i]
                TIJ = np.c_[self.ParticleType[n][i+1:], np.zeros_like(self.ParticleType[n][i+1:]) + self.ParticleType[n][i]]
                periodic = np.where(np.abs(RIJ / self.Boxlength) > 0.50, np.sign(RIJ), 0).astype(np.int)
                RIJ -= self.Boxlength * periodic * ppp    #remove PBC
                distance = np.sqrt(np.square(RIJ).sum(axis = 1))

                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 0] += Countvalue

                Countsum   = TIJ.sum(axis = 1)
                Countsub   = np.abs(TIJ[:, 0] - TIJ[:, 1])
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 1] += Countvalue #11
                Countvalue, BinEdge = np.histogram(distance[(Countsum == 4) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 2] += Countvalue #22
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 6], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 3] += Countvalue #33
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 4] += Countvalue #12
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 5], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 5] += Countvalue #23
                Countvalue, BinEdge = np.histogram(distance[(Countsum == 4) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 6] += Countvalue #13

        binleft  = BinEdge[:-1]   #real value of each bin edge, not index
        binright = BinEdge[1:]   #len(Countvalue) = len(BinEdge) - 1
        Nideal   = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
        grresults[:, 0]  = grresults[:, 0] * 2 / self.SnapshotNumber/ self.ParticleNumber / (Nideal * self.rhototal)
        grresults[:, 1]  = grresults[:, 1] * 2 / self.SnapshotNumber/ self.TypeNumber[0] / (Nideal * self.rhotype[0])
        grresults[:, 2]  = grresults[:, 2] * 2 / self.SnapshotNumber/ self.TypeNumber[1] / (Nideal * self.rhotype[1])
        grresults[:, 3]  = grresults[:, 3] * 2 / self.SnapshotNumber/ self.TypeNumber[2] / (Nideal * self.rhotype[2])
        grresults[:, 4]  = grresults[:, 4] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[1] /2.0
        grresults[:, 5]  = grresults[:, 5] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[1] / self.TypeNumber[2] /2.0
        grresults[:, 6]  = grresults[:, 6] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[2] /2.0

        binright = binright - 0.5 * rdelta #middle of each bin
        results  = np.column_stack((binright, grresults))
        names    = 'r  g(r)  g11(r)  g22(r)  g33(r)  g12(r)  g23(r)  g13(r)'
        np.savetxt(results_path + outputfile, results, fmt='%.6f', header = names, comments = '')
        print ('---------- Get g(r) results over ---------')
        return results

    def Quarternary(self, outputfile, ppp, rdelta = 0.01, results_path='../../analysis/gr/'):
        print ('--------- This is a Quarternary System ------')
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        grresults   = np.zeros((MAXBIN, 11))
        for n in range(self.SnapshotNumber):
            for i in range(self.ParticleNumber - 1):
                RIJ = self.Positions[n][i+1:] - self.Positions[n][i]
                TIJ = np.c_[self.ParticleType[n][i+1:], np.zeros_like(self.ParticleType[n][i+1:]) + self.ParticleType[n][i]]
                periodic = np.where(np.abs(RIJ / self.Boxlength) > 0.50, np.sign(RIJ), 0).astype(np.int)
                RIJ -= self.Boxlength * periodic * ppp    #remove PBC
                distance = np.sqrt(np.square(RIJ).sum(axis = 1))

                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 0] += Countvalue

                Countsum   = TIJ.sum(axis = 1)
                Countsub   = np.abs(TIJ[:, 0] - TIJ[:, 1])
                Countvalue, BinEdge = np.histogram(distance[Countsum == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 1] += Countvalue #11
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 2] += Countvalue #22
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 3] += Countvalue #33
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 8], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 4] += Countvalue #44
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 5] += Countvalue #12
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 6] += Countvalue #13
                Countvalue, BinEdge = np.histogram(distance[Countsub  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 7] += Countvalue #14
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 5) & (Countsub == 1)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 8] += Countvalue #23
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 9] += Countvalue #24
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 7], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:,10] += Countvalue #34

        binleft  = BinEdge[:-1]   #real value of each bin edge, not index
        binright = BinEdge[1:]   #len(Countvalue) = len(BinEdge) - 1
        Nideal   = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim)
        grresults[:, 0]  = grresults[:, 0] * 2 / self.SnapshotNumber/ self.ParticleNumber / (Nideal * self.rhototal)
        grresults[:, 1]  = grresults[:, 1] * 2 / self.SnapshotNumber/ self.TypeNumber[0] / (Nideal * self.rhotype[0])
        grresults[:, 2]  = grresults[:, 2] * 2 / self.SnapshotNumber/ self.TypeNumber[1] / (Nideal * self.rhotype[1])
        grresults[:, 3]  = grresults[:, 3] * 2 / self.SnapshotNumber/ self.TypeNumber[2] / (Nideal * self.rhotype[2])
        grresults[:, 4]  = grresults[:, 4] * 2 / self.SnapshotNumber/ self.TypeNumber[3] / (Nideal * self.rhotype[3])

        grresults[:, 5]  = grresults[:, 5] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[1] /2.0
        grresults[:, 6]  = grresults[:, 6] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[2] /2.0
        grresults[:, 7]  = grresults[:, 7] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[3] /2.0
        grresults[:, 8]  = grresults[:, 8] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[1] / self.TypeNumber[2] /2.0
        grresults[:, 9]  = grresults[:, 9] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[1] / self.TypeNumber[3] /2.0
        grresults[:,10]  = grresults[:,10] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[2] / self.TypeNumber[3] /2.0

        binright = binright - 0.5 * rdelta #middle of each bin
        results  = np.column_stack((binright, grresults))
        names    = 'r  g(r)  g11(r)  g22(r)  g33(r)  g44(r)  g12(r)  g13(r)  g14(r)  g23(r)  g24(r)  g34(r)'
        np.savetxt(results_path + outputfile, results, fmt='%.6f', header = names, comments = '')
        print ('---------- Get g(r) results over ---------')
        return results

    def Quinary(self, outputfile, ppp, rdelta = 0.01, results_path='../../analysis/gr/'):
        print ('--------- This is a Quinary System ------')
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        grresults   = np.zeros((MAXBIN, 16))
        for n in range(self.SnapshotNumber):
            for i in range(self.ParticleNumber - 1):
                RIJ = self.Positions[n][i+1:] - self.Positions[n][i]
                TIJ = np.c_[self.ParticleType[n][i+1:], np.zeros_like(self.ParticleType[n][i+1:]) + self.ParticleType[n][i]]
                periodic = np.where(np.abs(RIJ / self.Boxlength) > 0.50, np.sign(RIJ), 0).astype(np.int)
                RIJ -= self.Boxlength * periodic * ppp    #remove PBC
                distance = np.sqrt(np.square(RIJ).sum(axis = 1))

                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 0] += Countvalue

                Countsum   = TIJ.sum(axis = 1)
                Countsub   = np.abs(TIJ[:, 0] - TIJ[:, 1])
                Countvalue, BinEdge = np.histogram(distance[Countsum == 2], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 1] += Countvalue #11
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 2] += Countvalue #22
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 3] += Countvalue #33
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 8) & (Countsub == 0)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 4] += Countvalue #44
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 10], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 5] += Countvalue #55
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 3], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 6] += Countvalue #12
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 4) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 7] += Countvalue #13
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 5) & (Countsub == 3)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 8] += Countvalue #14
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 4)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:, 9] += Countvalue #15
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 5) & (Countsub == 1)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:,10] += Countvalue #23
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 6) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:,11] += Countvalue #24
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 7) & (Countsub == 3)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:,12] += Countvalue #25
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 7) & (Countsub == 1)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:,13] += Countvalue #34
                Countvalue, BinEdge = np.histogram(distance[(Countsum  == 8) & (Countsub == 2)], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:,14] += Countvalue #35
                Countvalue, BinEdge = np.histogram(distance[Countsum  == 9], bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults[:,15] += Countvalue #45

        binleft = BinEdge[:-1]   #real value of each bin edge, not index
        binright = BinEdge[1:]   #len(Countvalue) = len(BinEdge) - 1
        Nideal   = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim) 
        grresults[:, 0]  = grresults[:, 0] * 2 / self.SnapshotNumber/ self.ParticleNumber / (Nideal * self.rhototal)
        grresults[:, 1]  = grresults[:, 1] * 2 / self.SnapshotNumber/ self.TypeNumber[0] / (Nideal * self.rhotype[0])
        grresults[:, 2]  = grresults[:, 2] * 2 / self.SnapshotNumber/ self.TypeNumber[1] / (Nideal * self.rhotype[1])
        grresults[:, 3]  = grresults[:, 3] * 2 / self.SnapshotNumber/ self.TypeNumber[2] / (Nideal * self.rhotype[2])
        grresults[:, 4]  = grresults[:, 4] * 2 / self.SnapshotNumber/ self.TypeNumber[3] / (Nideal * self.rhotype[3])
        grresults[:, 5]  = grresults[:, 5] * 2 / self.SnapshotNumber/ self.TypeNumber[4] / (Nideal * self.rhotype[4])

        grresults[:, 6]  = grresults[:, 6] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[1] /2.0
        grresults[:, 7]  = grresults[:, 7] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[2] /2.0
        grresults[:, 8]  = grresults[:, 8] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[3] /2.0
        grresults[:, 9]  = grresults[:, 9] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[0] / self.TypeNumber[4] /2.0
        grresults[:,10]  = grresults[:,10] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[1] / self.TypeNumber[2] /2.0
        grresults[:,11]  = grresults[:,11] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[1] / self.TypeNumber[3] /2.0
        grresults[:,12]  = grresults[:,12] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[1] / self.TypeNumber[4] /2.0
        grresults[:,13]  = grresults[:,13] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[2] / self.TypeNumber[3] /2.0
        grresults[:,14]  = grresults[:,14] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[2] / self.TypeNumber[4] /2.0
        grresults[:,15]  = grresults[:,15] * 2 / self.SnapshotNumber/ Nideal * self.Volume /self.TypeNumber[3] / self.TypeNumber[4] /2.0

        binright = binright - 0.5 * rdelta #middle of each bin
        results  = np.column_stack((binright, grresults))
        names    = 'r  g(r)  g11(r)  g22(r)  g33(r)  g44(r)  g55(r)  g12(r)  g13(r)  g14(r)  g15(r)  g23(r)  g24(r)  g25(r)  g34(r)  g35(r)  g45(r)'
        np.savetxt(results_path + outputfile, results, fmt='%.6f', header = names, comments = '')
        print ('---------- Get g(r) results over ---------')
        return results

    def Senary(self, outputfile, ppp, rdelta = 0.01, results_path='../../analysis/gr/'):
        """ Only calculate the overall g(r) at this stage """
        print ('--------- This is a Senary System ------')
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        MAXBIN      = int(self.Boxlength.min() / 2.0 / rdelta)
        grresults   = np.zeros(MAXBIN)
        for n in range(self.SnapshotNumber):
            for i in range(self.ParticleNumber - 1):
                RIJ      = self.Positions[n][i+1:] - self.Positions[n][i]
                periodic = np.where(np.abs(RIJ / self.Boxlength) > 0.50, np.sign(RIJ), 0).astype(np.int)
                RIJ -= self.Boxlength * periodic * ppp    #remove PBC
                distance = np.sqrt(np.square(RIJ).sum(axis = 1))

                Countvalue, BinEdge = np.histogram(distance, bins = MAXBIN, range = (0, MAXBIN * rdelta))
                grresults += Countvalue

        binleft    = BinEdge[:-1]   #real value of each bin edge, not index
        binright   = BinEdge[1:]   #len(Countvalue) = len(BinEdge) - 1
        Nideal     = Nidealfac(self.ndim) * np.pi * (binright**self.ndim - binleft**self.ndim) 
        grresults  = grresults * 2 / self.ParticleNumber / self.SnapshotNumber / (Nideal * self.rhototal)

        binright = binright - 0.5 * rdelta #middle of each bin
        results  = np.column_stack((binright, grresults))
        names    = 'r  g(r)'
        np.savetxt(results_path + outputfile, results, fmt='%.6f', header = names, comments = '')
        print ('---------- Get g(r) results over ---------')
        return results




AlternativeAlgorithm = """
def Unary(self, outputfile, rdelta, ppp, results_path):        
    MAXBIN   = int(self.Boxlength.min() / 2.0 / rdelta)
    grresults = np.zeros(MAXBIN, dtype=np.int)

    for n in range(self.SnapshotNumber):
        for i in range(self.ParticleNumber - 1):
            RIJ = self.Positions[n][i+1:] - self.Positions[n][i]
            periodic = np.where(np.abs(RIJ / self.Boxlength) > 0.50, np.sign(RIJ), 0).astype(np.int)
            RIJ -= self.Boxlength * periodic * ppp    #remove PBC
            #RIJ -= np.dot(np.round(RIJ / self.Boxlength[np.newaxis, :]), np.diag(self.Boxlength)) * ppp #another way to remove PBC
            distance = np.sqrt(np.square(RIJ).sum(axis = 1))
            BIN = (distance / rdelta + 1).astype(np.int)
            BINCount = np.unique(BIN[BIN < MAXBIN], return_counts = True)  #This is very important due to numpy inherent characteristics
            grresults[BINCount[0]] += 2 * BINCount[1]

    f = open(results_path + outputfile, 'w')   #start to write results
    f.write('r  g(r)\n')                       #write file header
    for i in range(1, MAXBIN):                 #write data body
        Nideal = (4 / 3.0) * pi * (i**3 - (i-1)**3) * rdelta**3
        results = grresults[i] / self.SnapshotNumber / self.ParticleNumber / (Nideal * self.rhototal)
        f.write('{:.4f}  {:.6f}\n'.format((i-0.5)*rdelta, results))
    f.close()
"""
        