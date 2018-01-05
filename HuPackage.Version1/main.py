#coding = utf-8
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
#from dynamics import dynamics3d  #dynamics2d
# inputfile = './dump/dump.CuZr.2D.lammpstrj'
# dynoverall = 'overall.dat'
# dynpartial = 'partial.dat'
# dyn_result_path = './dynamics/'

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
