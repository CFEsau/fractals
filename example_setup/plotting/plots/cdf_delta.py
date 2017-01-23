#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import numpy as np
import sys
from sys import argv

snapnum, nmst, simnum, nran = argv[1:5] #snapnum passed from cdfplots.sh

nmst = int(nmst)
snapnum = int(snapnum)
nran=int(nran)

nedge = nmst-1

outpath = '../outputs/'
filepath = outpath+'runinv_k'+simnum+'/cluster_FoV5pc/'

# want differences between points at given edge number
# e.g. 3rd & 6th points for nmst=10
edgenum = [3,7]
# make list of delta x at edgenum[0] and edgenum[1]
# for all subsets of random stars (size 2 x nran):
dx = []
for i in range (2):#
    dx.append([])
    for j in range(0, nran):
        dx[i].append(0)
#print (dx)

#---------------------------
# CDF data for object stars
#---------------------------

objfile = 'MSTedges.dat'
objmacro = np.loadtxt(filepath+objfile)

# x values:
xobj=[0,0]

i=1
while i<nmst:
    if i==edgenum[0]:
        xobj[0]=objmacro[snapnum-1,i] #[,0] is snapnum, [,1] is l1, etc
        xobj[0]=float(xobj[0])
    if i==edgenum[1]:
        xobj[1]=objmacro[snapnum-1,i]
        xobj[1]=float(xobj[1])
    i=i+1
#print xobj[0:2]


#---------------------------
# CDF data for random stars
#---------------------------

#Number of plots for random MSTs:
j=1
while j<nran+1:
    
    filename = 'MSTedges_'+str(j)+'.dat'
    #data:
    ranmacro = np.loadtxt(filepath+filename)

    #x-values:
    xran=[0,0]
    i=1
    while i<nmst:
        if i==edgenum[0]:
            xran[0]=ranmacro[snapnum-1,i] #[,0] is snapnum, [,1] is l1, etc
            xran[0]=float(xran[0])
        if i==edgenum[1]:
            xran[1]=ranmacro[snapnum-1,i]
            xran[1]=float(xran[1])
        i=i+1
    #print filename,xran[0:2], xobj[0:2]
    
    dx[0][j-1]=float(xran[0]-xobj[0])
    dx[1][j-1]=float(xran[1]-xobj[1])
    #print 'dx:',dx[0][j-1],dx[1][j-1]
    
    j=j+1 #next file

# -----------------
# plot dx1 vs dx2
#------------------

#saveplot=''
plt.figure()
plt.xlabel("$\Delta$1")
plt.ylabel("$\Delta$2")
#find max dx value for axis limits:
axislim = max(np.max(dx),np.min(dx),key=abs)
plt.xlim(-axislim*1.05,axislim*1.05)
plt.ylim(-axislim*1.05,axislim*1.05)

#plot data:
plt.plot(dx[0][:],dx[1][:],linestyle="",marker="x")

#plt.show()
saveplot = (outpath+'plots/snap'+str(snapnum)+'_delta.png')
plt.savefig(saveplot, bbox_inches='tight')
print "Graph saved at " + saveplot
plt.close()
