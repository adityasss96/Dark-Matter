import pymultinest
import numpy as np
from scipy import sparse
from scipy.integrate import quad
import json
import sys
from numpy import log, exp, pi
import scipy.stats, scipy
import matplotlib.pyplot as plt

sigmazdata = None
SD=None
sigmarzdata = None
SDrz=None
densitydata=None
SDdensity=None
z=None
R=8
kdd=300

myfile=open('chains/1-sigma07')
f = myfile.read()
myfile.close()

#my_file=open('chains/1-sigma08')
#g = my_file.read().replace('\n',',')
#g = g.split(',')
#z = np.asarray(g)
#my_file.close()

den=float(f)


#def integrand2(z, R, kdd):
#      d = den*np.exp(-z/.9)*((-1*(2*267.65*z+(1500*z)/(np.sqrt(z**2+0.18**2))+(kdd*z)/(np.sqrt(z**2+2.5**2))))-(180.08*(z**1.44)*(1/R-(2/2.5))))
      #print(d)
#      return d

#def sigmaz():
#     x = np.linspace(0, z, 1000)
#     fx = np.exp(-x/.9)*((-1*(2*267.65*x+(1500*x)/(np.sqrt(x**2+0.18**2))+(kdd*x)/(np.sqrt(x**2+2.5**2))))-(180.08*(x**1.44)*(1/R-(2/2.5))))
#     area = np.sum(fx)*(z)/1000      
      
#     v = (area+1513.99892727)/(np.exp(-z/.9))     
#     return v**.5

def sigmarz(A):
      return A*(z**1.44)

def density():
      return den*(np.exp(-z/.9))


def prior(cube, ndim, nparams):
	cube[0] = (4*cube[0]*(10**4))**.5 
        
     
def loglikelihood(cube, ndim, nparams):
        A=cube[0]
        #sigmazmodel = sigmaz()
        sigmarzmodel = sigmarz(A)
        densitymodel = density()
        loglikelihood = (-0.5 *((sigmarzmodel - sigmarzdata) / SDrz)**2).sum() + (-0.5 *((densitymodel - densitydata) / SDdensity)**2).sum()
        
	return loglikelihood
	
parameters = ["A"]
n_params = len(parameters)

datafile = sys.argv[1]
sigmazdata = np.loadtxt(datafile)

datafile1 = sys.argv[2]
SD = np.loadtxt(datafile1)

datafile2 = sys.argv[3]
sigmarzdata = np.loadtxt(datafile2)

datafile3 = sys.argv[4]
SDrz = np.loadtxt(datafile3)

datafile4 = sys.argv[5]
densitydata = np.loadtxt(datafile4)

datafile5 = sys.argv[6]
SDdensity = np.loadtxt(datafile5)

datafile6 = sys.argv[7]
z = np.loadtxt(datafile6)


pymultinest.run(loglikelihood, prior, n_params, outputfiles_basename=datafile + '_1_', resume = False, verbose = False, n_live_points=1000)
json.dump(parameters, open(datafile + '_1_params.json', 'w')) 

a = pymultinest.Analyzer(outputfiles_basename=datafile + '_1_', n_params = n_params)

a_lnZ = a.get_stats()['global evidence']

print a_lnZ


