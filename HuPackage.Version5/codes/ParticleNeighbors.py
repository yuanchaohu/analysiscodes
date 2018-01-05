#!/usr/bin/python
# coding = utf-8
#This module is part of an analysis package

Authorinfo = """
             ------------------Name: Yuan-Chao Hu--------------
             --------------Email: ychu0213@gmail.com-----------
             ----------Web: https://yuanchaohu.github.io/------
             """

Docstr = """
             Reading partcles' Neighbor list and Voronoi polyhedron facearea 
             from the output of Voro++ Package analysis

             Voronoi tessellation can be carried out use the provided script 'Voronoi.sh'
             Voropp() is suitable for both data

         """

import numpy as np 

def Voropp(f, ParticleNumber):
    """
    Read Neighbor list data from the results of Voro++ Package
    &&&&&&&&&
    Read facearea list data from the results of Voro++ Package

    Read One Snapshot a time to save computer memory
    If you have multiple snapshots, you can import this function in a loop
    f = open(filename, 'f')
    The Voronoi analysis can be carried out use the provided shell secipt 'voronoi.sh'
    """

    header  = f.readline().split()  #header
    results = np.zeros((ParticleNumber, 50))

    for n in range(ParticleNumber):
        item = f.readline().split()
        results[int(item[0]) - 1, 0] = float(item[1])
        results[int(item[0]) - 1, 1:(int(item[1]) + 1)] = [float(j) - 1 for j in item[2:(int(item[1]) + 2)]]
        #Be attention to the '-1' after '=', all particle id has been reduced by 1
        #Please becareful when you are reading other data, like face area
        #you should first transform the data back before computation in other codes

    if 'neighborlist' in header:  #neighbor list should be integer
        results = results.astype(np.int)
    
    return results