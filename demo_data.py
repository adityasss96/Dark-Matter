import json
import sys
import numpy
from numpy import log, exp, pi
import scipy.stats, scipy
import matplotlib.pyplot as plt
import numpy as np
from random import uniform
from scipy.integrate import quad
from collections import defaultdict

R=8
kdd=300
no_data=10000
bins=20
d=defaultdict(list)
c=[]
w=[]
q=[]
f=defaultdict(list)
e=defaultdict(list)
g=defaultdict(list)
h=defaultdict(list)
m=defaultdict(list)
n=defaultdict(list)
j=[]
o=[]
y=[]
dd=defaultdict(list)
ee=defaultdict(list)
SD=defaultdict(list)
ts=defaultdict(list)
ss=defaultdict(list)
ds=defaultdict(list)

def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

def integrand1(z, R, kdd):
    d = np.exp(-z/.9)*(-(-1*(2*267.65*z+(1500*z)/(np.sqrt(z**2+0.18**2))+(kdd*z)/(np.sqrt(z**2+2.5**2))))+(180.08*(z**1.44)*(1/R-(2/2.5))))
    #print(d)
    return d

def e1(R, kdd):
    return quad(integrand1, 0, np.inf, args=(R, kdd))[0]

def integrand2(z, R, kdd):
      d = np.exp(-z/.9)*((-1*(2*267.65*z+(1500*z)/(np.sqrt(z**2+0.18**2))+(kdd*z)/(np.sqrt(z**2+2.5**2))))-(180.08*(z**1.44)*(1/R-(2/2.5))))
      #print(d)
      return d

def sigmaz(z):
      def e2(R, kdd):
       return quad(integrand2, 0, z, args=(R, kdd))[0]

      d = (e2(R, kdd)+e1(R, kdd))/(np.exp(-z/.9))     
      return d
def sigmarz(z):
      return 180.08*(z**1.44)


for i in range(no_data):
        z = np.random.exponential(scale=0.9, size=None)
        q.append(z)	
x=np.hstack(q)        

for l in range(no_data):
 
 b = histedges_equalN(x, 20)
 
 for p in range(20):
   if b[p]<=x[l]<b[p+1]:
     a=x[l]
     sigmazdata = np.random.normal(loc=a,scale=sigmaz(a)**.5)   
     f[p]=sigmazdata
     e[p].append(f[p])
     g[p]=x[l]
     h[p].append(g[p])
     sigmarzdata = np.random.normal(loc=a,scale=sigmarz(a))
     m[p]=sigmarzdata
     n[p].append(m[p])
  
   p=p+1
 l=l+1

data=np.concatenate((np.hstack(e[0]),np.hstack(e[1]),np.hstack(e[2]),np.hstack(e[3]),np.hstack(e[4]),np.hstack(e[5]),np.hstack(e[6]),np.hstack(e[7]),np.hstack(e[8]),np.hstack(e[9]),np.hstack(e[10]),np.hstack(e[11]),np.hstack(e[12]),np.hstack(e[13]),np.hstack(e[14]),np.hstack(e[15]),np.hstack(e[16]),np.hstack(e[17]),np.hstack(e[18]),np.hstack(e[19])),axis=0)

numpy.savetxt("chains/1-sigma%02d" % (1), data )

data2=np.concatenate((np.hstack(n[0]),np.hstack(n[1]),np.hstack(n[2]),np.hstack(n[3]),np.hstack(n[4]),np.hstack(n[5]),np.hstack(n[6]),np.hstack(n[7]),np.hstack(n[8]),np.hstack(n[9]),np.hstack(n[10]),np.hstack(n[11]),np.hstack(n[12]),np.hstack(n[13]),np.hstack(n[14]),np.hstack(n[15]),np.hstack(n[16]),np.hstack(n[17]),np.hstack(n[18]),np.hstack(n[19])),axis=0)

numpy.savetxt("chains/1-sigma%02d" % (3), data2 )

data=np.concatenate((np.hstack(h[0]),np.hstack(h[1]),np.hstack(h[2]),np.hstack(h[3]),np.hstack(h[4]),np.hstack(h[5]),np.hstack(h[6]),np.hstack(h[7]),np.hstack(h[8]),np.hstack(h[9]),np.hstack(h[10]),np.hstack(h[11]),np.hstack(h[12]),np.hstack(h[13]),np.hstack(h[14]),np.hstack(h[15]),np.hstack(h[16]),np.hstack(h[17]),np.hstack(h[18]),np.hstack(h[19])),axis=0)

numpy.savetxt("chains/1-sigma%02d" % (8), data )

for i in range(19):
   for p in range(no_data/20): 
     std=np.std(np.hstack(e[i]))
     c.append(std)
     p+=1
   i+=1
for i in range((no_data/20)-1):
   std=np.std(np.hstack(e[19]))
   c.append(std)
   i+=1
data1=np.hstack(c) 

numpy.savetxt("chains/1-sigma%02d" % (2), data1 )

for i in range(19):
   for p in range(no_data/20): 
     std=np.std(np.hstack(n[i]))
     w.append(std)
     p+=1
   i+=1
for i in range((no_data/20)-1):
   std=np.std(np.hstack(n[19]))
   w.append(std)
   i+=1
data3=np.hstack(w) 
numpy.savetxt("chains/1-sigma%02d" % (4), data3 )


for j in range(20):
  cc=[]
  for i in range(no_data):
    aa = np.random.exponential(scale=0.9, size=None)
    i = i+1
    cc.append(aa)
    ts[j].append(aa) 
  s=np.hstack(cc)
  ss[j]=np.hstack(ts[j])
  
  for k in range(bins): 
   width_bins = histedges_equalN(s, 20)[k+1]-histedges_equalN(s, 20)[k]
   dd[k] = (no_data/bins)/(3.14*width_bins)    
   ee[k].append(dd[k])
   k = k+1
   
  j = j+1  
  
for l in range(bins):
 SD[l]=np.std(np.hstack(ee[l]))
 l=l+1

for i in range(20):
  for p in range(20):
    width_bins = histedges_equalN(ss[i], 20)[p+1]-histedges_equalN(ss[i], 20)[p]
    d[p] = (no_data/20)/(3.14*width_bins)
    ds[p].append(d[p])
    p+=1


for i in range(19):
   for p in range(no_data/20): 
     y.append(np.mean(np.hstack(ds[i])))
     p+=1
   i+=1
for i in range((no_data/20)-1):
   y.append(np.mean(np.hstack(ds[19])))
   i+=1
data4=np.hstack(y) 
numpy.savetxt("chains/1-sigma%02d" % (5), data4 )


for i in range(19):
   for p in range(no_data/20): 
     std=SD[i]
     o.append(std)
     p+=1
   i+=1
for i in range((no_data/20)-1):
   std=SD[19]
   o.append(std)
   i+=1
data5=np.hstack(o) 
numpy.savetxt("chains/1-sigma%02d" % (6), data5 )

np.savetxt("chains/1-sigma%02d" % (7), (np.mean(np.hstack(ee[0]))+2*np.std(np.hstack(ee[0]))).reshape(1,))


for i in range():
    x = np.linspace(0, z, 1000)
    fx = np.exp(-x/.9)*((-1*(2*267.65*x+(1500*x)/(np.sqrt(x**2+0.18**2))+(kdd*x)/(np.sqrt(x**2+2.5**2))))-(180.08*(x**1.44)*(1/R-(2/2.5))))
    area = np.sum(fx)*(z)/1000 


