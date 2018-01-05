##!/usr/bin/python
# coding = utf-8
#This module is part of an analysis package

Authorinfo = """
             ------------------Name: Yuan-Chao Hu--------------
             --------------Email: ychu0213@gmail.com-----------
             ----------Web: https://yuanchaohu.github.io/------
             """

Docstr = """
         This module is responsible for performing Voronoi tessellation 
         By using the Voro++ Package (http://math.lbl.gov/voro++/about.html)
         """

import numpy  as np 
import pandas as pd
import os, subprocess
from dump import readdump


def get_input(inputfile, ndim, radii):
    """ Design input file for Voro++ by considering particle radii 
        radii must be a dict like {1 : 1.28, 2 : 1.60}
        if you do not want to consider radii, set the radii the same
    """

    d = readdump(inputfile, ndim)
    d.read_onefile()
    results = []
    for n in range(d.SnapshotNumber):
        ParticleRadii = np.array(pd.Series(d.ParticleType[n]).map(radii))
        PositionRadii = np.column_stack((d.Positions[n], ParticleRadii))
        voroinput     = np.column_stack((np.arange(d.ParticleNumber[n]) + 1, PositionRadii))
        results.append(voroinput)

    return (results, d.Boxbounds)

def cal_voro(inputfile, ndim, radii, ppp = '-p', results_path = '../../analysis/voro/'):
    """ Radical Voronoi Tessellation using voro++ 

        radii must be a dict like {1 : 1.28, 2 : 1.60}
        if you do not want to consider radii, set the radii the same
        There are two methods in choosing box boundaries
        One is from the inherent snapshot
        The other from the minimum and maximum of particle coordinates
        The results are influenced by this choice
        Set ppp as '-p' for periodic boundary conditions at all direction
        Set ppp for each direction from '-px -py -pz' 
    """   
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    basename  = os.path.splitext(os.path.basename(inputfile))[0]
    fneighbor = open(results_path + basename + '.neighbor.dat', 'w')
    ffacearea = open(results_path + basename + '.facearea.dat', 'w')
    findex    = open(results_path + basename + '.voroindex.dat', 'w') 
    foverall  = open(results_path + basename + '.overall.dat', 'w')

    position, bounds = get_input(inputfile, ndim, radii)
    for n in range(len(position)):
        fileformat = '%d ' + '%.6f ' * ndim + '%.2f'
        np.savetxt('dumpused', position[n], fmt = fileformat)
        
        #use box boundaries from snapshot
        Boxbounds = bounds[n].ravel()
        #use box boundaries from particle coordinates 
        # boundsmin = position[n][:, 1: ndim + 1].min(axis = 0) - 0.1
        # boundsmax = position[n][:, 1: ndim + 1].max(axis = 0) + 0.1
        # Boxbounds = (np.column_stack((boundsmin, boundsmax))).ravel()

        cmdline = 'voro++ ' + ppp + ' -r -c "%i %s %v %F @%i %A @%i %s %n @%i %s %f" '\
                  + ('%f %f ' * ndim % tuple(Boxbounds)) + 'dumpused'
        if n == 0: print (cmdline)
        subprocess.run(cmdline, shell = True)

        fneighbor.write('id   cn   neighborlist\n')
        ffacearea.write('id   cn   facearealist\n')
        findex.write('id   voro_index   0_to_7_faces\n')
        foverall.write('id   cn   volume   facearea\n')
        f = open('dumpused.vol', 'r')
        for i in range(len(position[n][:, 0])):
            item = f.readline().split('@')
            foverall.write(item[0]  + '\n')
            findex.write(item[1] + '\n')
            fneighbor.write(item[2] + '\n')
            ffacearea.write(item[3])
        f.close()

    os.remove('dumpused')      #delete temporary files
    os.remove('dumpused.vol')
    fneighbor.close()
    ffacearea.close()
    foverall.close()
    findex.close()
    print ('---------- Voronoi Analysis Done ------------')