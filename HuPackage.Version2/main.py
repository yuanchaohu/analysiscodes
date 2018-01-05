#coding = utf8
import sys
import numpy as np 
sys.path.append('./codes/')

#--------------------Read Snapshots-------------
#from dump import readdump
#inputfile = './dumpfile'
#readdump(inputfile)  #get all the informtion in dump files (one file)
#x, xs, xu coordinates are accepted. 
#This part is not necesay since dump has been imported in follow modules


#--------------------Calculate Pair Correlation Functions-------------
#from paircorrelationfunctions import gr3d  #gr2d
#inputfile = './dumpfile'
#gr3d(inputfile).getresults(outputfile = 'gr.dat', results_path = './gr/')

#gr3d(inputfile).getresults(outputfile, rdelta = 0.01, ppp = [1,1,1], results_path = '../../analysis/gr/')
#gr2d(inputfile).getresults(outputfile, rdelta = 0.01, ppp = [1,1], results_path = '../../analysis/gr/')
#outputfile: file name of output results without path
#rdelta: bin value for gr, default = 0.01 (2D & 3D)
#ppp: set periodic boundary conditions, 1 for yes 0 for no. default = [1,1,1] (3D) and [1,1] (2D)
#results_path: set path for outcomes, default = '../../analysis/gr/' (2D & 3D)
#Individual functions Unary(), Binary(), Ternary(), Quarternary(), Qinary(), Senary() can also be imported the same as getresults()
#A numpy array of output results is returned in individual functions other than getresults()
#gr3d(sys.argv[1]).getresults(outputfile = sys.argv[2], rdelta = sys.argv[3], ppp = sys.argv[4], results_path = sys.argv[5]) for shell


#--------------------Calculate Structure Factors-------------
#from structurefactors import sq3d  #sq2d
#inputfile = './dumpfile'
#sqoutputfile = 'testUnary.dat'
#sq_results_path = './Sq/'
#sq3d(inputfile).getresults(outputfile = sqoutputfile, results_path = sq_results_path)

#sq3d(inputfile).getresults(outputfile, results_path = '../../analysis/sq/')
#sq2d(inputfile).getresults(outputfile, results_path = '../../analysis/sq/')
#outputfile: file name of output results without path
#results_path: set path for outcomes, default = '../../analysis/sq/' (2D & 3D)
#Individual functions Unary(), Binary(), Ternary(), Quarternary(), Qinary(), Senary() can also be imported the same as getresults()
#A numpy array of output results is returned in the individual functions other than getresults()

#from structurefactors import wavevector3d  #wavevector2d
#wavevector2d/3d(Numofq = 500) 
#function setting wave vector for structure factors
#return wave vector numpy array [d, a, b, c] d = a**2 + b**2 + c**2 (3D)
#return wave vector numpy array [d, a, b] d = a**2 + b**2 (2D)
#Numofq: the number of wavenumber considered, default = 500


# #--------------------Calculate Dynamical Properties-------------
# from dynamics import dynamics3d  #dynamics2d
# inputfile = './dump/CuZr.xu.lammpstrj'
# dyn_result_path = './dynamics/'
# dynamics3d(inputfile).slowS4(outputfile = 'testS4slow.Binary.dat', X4time = 0.5, dt = 0.002, a = 1.0, results_path = dyn_result_path)
# dynamics3d(inputfile).fastS4(outputfile = 'testS4fast.Binary.dat', a = 1.0, dt = 0.002, X4timeset = 0.5, results_path = dyn_result_path)


#get the first peaks of structure factors first for self-intermediate scattering functions
#sqvalue = np.loadtxt(sq_results_path + sqoutputfile, skiprows = 1)
#qmaxvalue = [sqvalue[sqvalue[:, i].argmax(), 0] for i in range(1, len(sqvalue[0]))]
#dynamics3d(inputfile).total(outputfile = dynoverall, qmax = qmaxvalue[0], a = 1.0, dt = 0.002, results_path = '../../analysis/dynamics'):
#dynamics3d(inputfile).partial(outputfile = dynpartial, qmax = qmaxvalue[1:], a = 1.0, dt = 0.002, results_path = '../../analysis/dynamics'):
#qmax should be a list for partial(), such as [0, 0, 0] for ternary system. A list is returned containing results of different types in squence
#getX4value = np.loadtxt(dyn_result_path + dynoverall, skiprows = 1)
#getX4peak = getX4value[getX4value[:, 4].argmax(), 0] #just the peak time value not the interval
#dynamics3d(inputfile).slowS4(outputfile, X4time = getX4peak, dt = 0.002, a = 1.0, results_path = '../../analysis/dynamics'):
#dynamics3d(inputfile).fastS4(outputfile, a = 1.0, dt = 0.002, X4timeset = 0, results_path = '../../analysis/dynamics'):
#X4timeset can be set as in slowS4(), but can also be 0 to use the inner calculated value in fastS4(). default = 0

#get the first peaks of structure factors first for self-intermediate scattering functions
#sqvalue = np.loadtxt(sq_results_path + sqoutputfile, skiprows = 1)
#qmaxvalue = [sqvalue[sqvalue[:, i].argmax(), 0] for i in range(1, len(sqvalue[0]))]
#dynamics2d(inputfile).total(outputfile = dynoverall, qmax = qmaxvalue[0], a = 1.0, dt = 0.002, results_path = '../../analysis/dynamics'):
#dynamics2d(inputfile).partial(outputfile = dynpartial, qmax = qmaxvalue[1:], a = 1.0, dt = 0.002, results_path = '../../analysis/dynamics'):
#qmax should be a list for partial(), such as [0, 0, 0] for ternary system. A list is returned containing results of different types in squence
#getX4value = np.loadtxt(dyn_result_path + dynoverall, skiprows = 1)
#getX4peak = getX4value[getX4value[:, 4].argmax(), 0] #just the peak time value not the interval
#dynamics2d(inputfile).slowS4(outputfile, X4time = getX4peak, dt = 0.002, a = 1.0, results_path = '../../analysis/dynamics'):
#dynamics2d(inputfile).fastS4(outputfile, a = 1.0, dt = 0.002, X4timeset = 0, results_path = '../../analysis/dynamics'):
#X4timeset can be set as in slowS4(), but can also be 0 to use the inner calculated value in fastS4(). default = 0

#outputfile: file name of output results without path
#Function total() only calculates the overall particle dynamics neglecting particle types
#Use Function partial() to calculate particle type related dynamics
#Use Function slowS4() to calculate four-point dynamic structure factor of slow atoms
#Use Function fastS4() to calculate four-point dynamic structure factor of fast atoms
#a: cutoff in the overlap function Q(t). default = 1.0
#dt: timestep in MD simulations. default = 0.002

# #--------------------Calculate Cage Relative Dynamical Properties-------------
#from dynamics import dynamics3d  #dynamics2d
#inputfile = './dump/ZrCuAl.xu.lammpstrj'
#dyn_result_path = './dynamics/'
#dynamics3d(inputfile).total(outputfile = 'total.ZrCuAl.dat', qmax = 2.5, a = 1.0, dt = 0.002, results_path = dyn_result_path)
#dynamics3d(inputfile).partial(outputfile = 'partial.ZrCuAl.dat', qmax = [2.78, 2.56, 2.5], a = 1.0, dt = 0.002, results_path = dyn_result_path)
#dynamics3d(inputfile).slowS4(outputfile = 'S4slow.ZrCuAl.dat', X4time = 0.7, dt = 0.002, a = 1.0, results_path = dyn_result_path)
#dynamics3d(inputfile).fastS4(outputfile = 'S4fast.ZrCuAl.dat', a = 1.0, dt = 0.002, X4timeset = 0.7, results_path = dyn_result_path)
#the synax of this module is the same as above of the Dynamical Properties
#but the neighbor list file should also be provided to calculate the cage relative displacements

# #--------------------Calculate BOO in 3D -------------
# from BondOOrder import BOO3D
# for j in ['hcp', 'fcc']: #['sc', 'bcc', 'fcc', 'hcp']
#     Neighborfile = './neighbor/' + j + '.neighborlist.dat'
#     faceareafile = './neighbor/' + j + '.facearealist.dat'
#     dumpfile     = './dump/dump.' + j + '.lammpstrj'
#     filepath = './boo/' + j +'/'
#     print (filepath)
#     boo = BOO3D(dumpfile, Neighborfile, faceareafile)
#     for lll in range (10, 11): #(4, 11, 2):
#         print (lll)
#         #boo.qlQl(l = lll, ppp = [1,1,1], AreaR = 0, outputql = 'sq' + str(lll) +'.dat', outputQl = 'bQ' + str(lll) +'.dat', results_path = filepath)
#         boo.smallwcap(l = lll, ppp = [1,1,1], AreaR = 0, outputw = 'sw' + str(lll) +'.dat', outputwcap = 'swcap' + str(lll) +'.dat', results_path = filepath)

# boo = BOO3D(dumpfile = './dump/CuZr.x.lammpstrj', Neighborfile = './neighbor/voro.CuZr.x.neighbor.dat', faceareafile = './neighbor/voro.CuZr.x.facearea.dat')
# boo.qlQl(l = lll, ppp = [1,1,1], AreaR = 0, outputql = 'sq' + str(lll) +'.dat', outputQl = 'bQ' + str(lll) +'.dat', results_path = filepath)
# boo.sijsmallql(l = lll, ppp = [1,1,1], AreaR = 0, c = 0.7, outputql = 'sijsq' + str(lll) +'.dat', outputsij = 'sijvaluesq' + str(lll) +'.dat', results_path = filepath)
# boo.sijlargeQl(l = lll, ppp = [1,1,1], AreaR = 0, c = 0.7, outputQl = 'sijbQ' + str(lll) +'.dat', outputsij = 'sijvaluebQ' + str(lll) +'.dat', results_path = filepath)
# boo.GllargeQ(l = lll, ppp = [1,1,1], rdelta = 0.01, AreaR = 0, outputgl = 'GlbQ' + str(lll) +'.dat', results_path = filepath)
# boo.Glsmallq(l = lll, ppp = [1,1,1], rdelta = 0.01, AreaR = 0, outputgl = 'Glsq' + str(lll) +'.dat', results_path = filepath)
# boo.smallwcap(l = lll, ppp = [1,1,1], AreaR = 0, outputw = 'sw' + str(lll) +'.dat', outputwcap = 'swcap' + str(lll) +'.dat', results_path = '../../analysis/BOO')
# boo.largeWcap(l = lll, ppp = [1,1,1], AreaR = 0, outputW = 'bW' + str(lll) +'.dat', outputWcap = 'bWcap' + str(lll) +'.dat', results_path = '../../analysis/BOO')
# boo.timecorr(l = lll, ppp = [1,1,1], AreaR = 0, dt = 0.002, outputfile = 'timecorr' + str(lll) +'.dat', results_path = filepath)


# boo = BOO3D(dumpfile = , Neighborfile = , faceareafile = )
# boo.qlQl(l = , ppp = [1,1,1], AreaR = 0, outputql = '', outputQl = '', results_path = '../../analysis/BOO')
# #Give names to outputql and outputQl to store the data or no data will be dumped
# #Both results will be returned in a tuple (ql, Ql)
# boo.sijsmallql(l = , ppp = [1,1,1], AreaR = 0, c = 0.7, outputql = '', outputsij = '', results_path = '../../analysis/BOO')
# boo.sijlargeQl(l = , ppp = [1,1,1], AreaR = 0, c = 0.7, outputQl = '', outputsij = '', results_path = '../../analysis/BOO')
# #c is a cutoff demonstrating whether a bond is crystalline or not
# #Give names to outputql and outputsij to store the results
# #only individual sij results will be returned as a numpy array
# boo.GllargeQ(l = , ppp = [1,1,1], rdelta = 0.01, AreaR = 0, outputgl = 'GllargeQ.dat', results_path = '../../analysis/BOO')
# boo.Glsmallq(l = , ppp = [1,1,1], rdelta = 0.01, AreaR = 0, outputgl = 'Glsmallq.dat', results_path = '../../analysis/BOO')
# #rdelta is the bin size in calculating g(r) and Gl(r)
# #results of g(r) and Gl(r) are returned as a numpy array
# boo.BoowWcap(l = , ppp = [1,1,1], AreaR = 0, outputw = '', outputW = '', outputwcap = '', outputWcap = '', results_path = '../../analysis/BOO')
#Give names to outputw, outputW, outputwcap, outputWcap to store the data, default is none
#results of w and W including normalized ones are returned in a tuple (w, W, w_cap, W_cap)
#timecorr(l, ppp = [1,1,1], AreaR = 0, dt = 0.002, outputfile = 'timecorr.dat', results_path = '../../analysis/BOO'):
#AreaR = 0 indicates calculate traditional ql and Ql.  default = 0
#AreaR = 1 indicates calculate voronoi polyhedron face area weighted ql and Ql
#ppp: set periodic boundary conditions, 1 for yes 0 for no. default = [1,1,1] (3D)

# --------------------Rewrite Neighbor list from LAMMPS Voronoi -------------
# from ParseList import readall
# for j in ['sc', 'bcc', 'fcc', 'hcp']:
#     fnfile = './neighbor/' + j + '.neighborlist.dat'
#     fffile = './neighbor/' + j + '.facearealist.dat'
#     fread  = './dump/dump.' + j + '.neighbors.dat'
#     if j is 'sc': ParticleNumber = 1000
#     if j is 'bcc': ParticleNumber = 2000
#     if j is 'fcc': ParticleNumber = 4000
#     if j is 'hcp': ParticleNumber = 4000
#     readall(fnfile, fffile, fread, ParticleNumber, Snapshotnumber = 2)

#--------------------Calculate Binder cumulant -------------
# from Binder import BCumulant
# BCumulant('./dump/CuZr.x.lammpstrj', 3).BCall(outputfile = 'binderCuZr.dat', results_path = './dynamics/')

#--------------------Calculate Bond Orietational Order at 2D -------------
from Order2D import BOO2D

dumpfile = './dump/dump.CuZr.2D.lammpstrj'
Neighborfile = './Order2D/neighborlist-gr-0K-1.dat'
BOO2D(dumpfile, Neighborfile).tavephi(outputphi = 'phi6.dat', outputavephi = 'ave_phi6.dat', l = 6, ppp = [1, 1], avet = 64, dt = 0.002, results_path = './Order2D/')
BOO2D(dumpfile, Neighborfile).spatialcorr(outputfile = 'g6.dat', l = 6, ppp = [1, 1], rdelta = 0.01, results_path = './Order2D/')
BOO2D(dumpfile, Neighborfile).timecorr(outputfile = 'time6.dat', l = 6, ppp = [1, 1], dt = 0.002, results_path = './Order2D/')

# BOO2D(dumpfile, Neighborfile).tavephi(outputphi, outputavephi, avet, l = 6, ppp = [1, 1], dt = 0.002, results_path = '../../analysis/order2d/')
# BOO2D(dumpfile, Neighborfile).spatialcorr(outputfile, l = 6, ppp = [1, 1], rdelta = 0.01, results_path = '../../analysis/order2d/')
# BOO2D(dumpfile, Neighborfile).timecorr(outputfile, l = 6, ppp = [1, 1], dt = 0.002, results_path = '../../analysis/order2d/')
