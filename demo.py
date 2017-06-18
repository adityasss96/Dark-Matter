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
area=None
z=None
#R=8
#kdd=300

myfile=open('chains/1-sigma07')
f = myfile.read()
myfile.close()

den=float(f)


#1513.99892727
def sigmaz(A,h,n):
     
     pp=[]
     
     def integrand2(z):
        dd = np.exp(-z/h)*(-1*(-1*(2*267.65*z+(1500*z)/(np.sqrt(z**2+0.18**2))+(300*z)/(np.sqrt(z**2+2.5**2))))+(A*(z**n)*(1/8-(2/2.5))))
       
        return dd
     uptoinf = scipy.integrate.quad(integrand2, 0, np.inf)[0]
   
     def integrand1(x):
       
        d = np.exp(-x/h)*((-1*(2*267.65*x+(1500*x)/(np.sqrt(x**2+0.18**2))+(300*x)/(np.sqrt(x**2+2.5**2))))-(A*(x**n)*(1/8-(2/2.5))))
        return d 
        

     for i in range(len(z)): 
      
      x=z[i]

      inte = scipy.integrate.quad(integrand1, 0, z[i])[0]
      
      """
      x = np.linspace(0, z[i], 100000000)
      
      y= np.exp(-x/h)*((-1*(2*267.65*x+(1500*x)/(np.sqrt(x**2+0.18**2))+(kdd*x)/(np.sqrt(x**2+2.5**2))))-(A*(x**n)*(1/R-(2/2.5))))
     
      inte = np.trapz(y,x) 
      """
      pp.append(inte)
      
     #print np.hstack(pp)
     v = (np.hstack(pp)+uptoinf)/(np.exp(-z/h))
     #print z
     #print np.hstack(pp)  
     #print v
     #print v**.5
     #print
     return v**.5

def sigmarz(A,n):
      return A*(z**n)

def density(h):
      return den*(np.exp(-z/h))


def prior(cube, ndim, nparams):
	cube[0] = (4*cube[0]*(10**4))**.5 
        cube[1] = (((1.9**2)-1)*cube[1]+1)**.5
        cube[2] = (((1.4**2)-(.4**2))*cube[2]+(.4**2))**.5     

def loglikelihood(cube, ndim, nparams):
        
        A=cube[0]
        n=cube[1]
        h=cube[2]
        sigmazmodel = sigmaz(A,h,n)
        sigmarzmodel = sigmarz(A,n)
        densitymodel = density(h)
        loglikelihood = (-0.5 *((sigmazmodel - sigmazdata) / SD)**2).sum() + (-0.5 *((sigmarzmodel - sigmarzdata) / SDrz)**2).sum() + (-0.5 *((densitymodel - densitydata) / SDdensity)**2).sum()
        
	return loglikelihood
	
parameters = ["A", "n", "h"]
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

datafile7 = sys.argv[8]
area = np.loadtxt(datafile7)


pymultinest.run(loglikelihood, prior, n_params, outputfiles_basename=datafile + '_1_', resume = False, verbose = False, n_live_points=400)
json.dump(parameters, open(datafile + '_1_params.json', 'w')) 

a = pymultinest.Analyzer(outputfiles_basename=datafile + '_1_', n_params = n_params)

a_lnZ = a.get_stats()['global evidence']

print a_lnZ


