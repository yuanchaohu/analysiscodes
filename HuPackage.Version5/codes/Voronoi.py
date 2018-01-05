##!/usr/bin/python
# coding = utf-8
#This module is part of an analysis package

Authorinfo = """
             ------------------Name: Yuan-Chao Hu--------------
             --------------Email: ychu0213@gmail.com-----------
             ----------Web: https://yuanchaohu.github.io/------
             """

Docstr = """
         This module is responsible for performing (radical) Voronoi tessellation 
         By using the Voro++ Package (http://math.lbl.gov/voro++/about.html)

         Use cal_voro() for periodic boundary conditions
         Use voronowalls() for non-periodic boundary conditions because there are artifacial walls
         Use indicehis() to statistics the Voronoi Index
         """

import numpy  as np 
import pandas as pd
import os, subprocess, re 
from   dump import readdump


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
    """ Radical Voronoi Tessellation using voro++ Originally

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
    findex.write('id   voro_index   0_to_7_faces\n')
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

def voronowalls(inputfile, ndim, radii, ppp, results_path = '../../analysis/voro/'):
    """ Radical Voronoi Tessellation using voro++ 
        Output Results by Removing Artifacial Walls

        radii must be a dict like {1 : 1.28, 2 : 1.60}
        if you do not want to consider radii, set the radii the same
        There are two methods in choosing box boundaries
        One is from the inherent snapshot
        The other from the minimum and maximum of particle coordinates
        The results are influenced by this choice
        Set ppp for each direction from '-px -py -pz' 
    """   
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    basename  = os.path.splitext(os.path.basename(inputfile))[0]
    fneighbor = open(results_path + basename + '.neighbor.dat', 'w')
    ffacearea = open(results_path + basename + '.facearea.dat', 'w')
    findex    = open(results_path + basename + '.voroindex.dat', 'w') 
    findex.write('id   voro_index   0_to_7_faces\n')
    foverall  = open(results_path + basename + '.overall.dat', 'w')
    np.set_printoptions(threshold = np.inf, linewidth = np.inf)

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
        foverall.write('id   cn   volume   facearea\n')
        f = open('dumpused.vol', 'r')
        for i in range(len(position[n][:, 0])):
            item        = f.readline().split('@')

            medium      = [int(j)   for j in item[2].split()]
            mneighbor   = np.array(medium, dtype = np.int)
            neighbor    = mneighbor[mneighbor > 0]
            neighbor[1] = len(neighbor[2:])

            medium      = [float(j) for j in item[3].split()]
            facearea    = np.array(medium)
            facearea    = facearea[mneighbor > 0]
            facearea[1] = neighbor[1]

            medium      = [float(j) for j in item[0].split()]
            overall     = np.array(medium)
            overall[1]  = neighbor[1]
            overall[3]  = facearea[2:].sum()

            #-----write Overall results-----
            np.savetxt('temp', overall[np.newaxis, :], fmt = '%d %d %.6f %.6f')
            with open('temp') as temp:
                foverall.write(temp.read())
            #-----write voronoi index-------
            findex.write(item[1] + '\n')
            #-----write facearea list-------
            np.savetxt('temp', facearea[np.newaxis, :], fmt = '%d ' * 2 + '%.6f ' * neighbor[1])
            with open('temp') as temp:
                ffacearea.write(temp.read())
            #-----write neighbor list-------
            fneighbor.write(re.sub('[\[\]]', ' ', np.array2string(neighbor) + '\n'))
        f.close()

    os.remove('dumpused')      #delete temporary files
    os.remove('dumpused.vol')
    os.remove('temp')
    fneighbor.close()
    ffacearea.close()
    foverall.close()
    findex.close()
    print ('---------- Voronoi Analysis Done ------------')

def indicehis(inputfile, outputfile, SnapshotNumber, ParticleNumber, results_path = '../../analysis/voro/'):
    """ Statistics of Voronoi Indices to Give the Frequency """
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    f = open(inputfile,  'r')
    f.readline()
    totaldata = SnapshotNumber * ParticleNumber
    medium = np.zeros((totaldata, 12), dtype = np.int)
    for n in range(totaldata):
        item = f.readline().split()[1:]
        medium[n, :len(item)] = item #[int(j) for j in item]

    medium          = medium[:, 2:8]
    indices, counts = np.unique(medium, return_counts = True, axis = 0)
    sort_indices    = np.argsort(counts)[::-1]
    freq_indices    = counts / ParticleNumber / SnapshotNumber
    results         = np.column_stack((indices[sort_indices], freq_indices[sort_indices]))
    fformat         = '%d ' * 6 + '%.6f '
    names           = 'Voronoi Indices,  Frequency'
    np.savetxt(results_path + outputfile, results[:50], fmt = fformat, header = names, comments = '')
    f.close()
    return results

