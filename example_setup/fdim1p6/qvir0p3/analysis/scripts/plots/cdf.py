#!/usr/bin/env python


import os
import matplotlib.pyplot as plt
import numpy as np
import sys


snapnum=70 #plot the cdf for this snapshot
nmst=10
nedge=nmst-1

#print ("   Plotting CMD...", sep="")

my_dpi=96
saveplot = ''

plt.figure()
plt.xlabel("Edge")
plt.ylabel("Cumulative edge length")
plt.title("CDF for MST edge lengths")
plt.xlim(xmax = nedge, xmin = 0)
#plt.ylim(ymax = 1, ymin = 0)
filepath = '../../outputs/runinv_k01/cluster_all/MSTedges_xy.dat'

#data:
macro = np.loadtxt(filepath)


#x-values:

iedge=np.arange(1,nmst,1) #this gives integer numbers from 1 to nedge


#y-values:
#want to plot cumulative edgelengths 1:nmst-1 for given snapshot

#First y value is just smallest edge length from 'macro' array:]
edgelength=[]
edgelength.append(macro[snapnum-1,1])

#Additional edge lengths are cumulative:
i=1
while i<nedge:
    newedge=macro[snapnum-1,i+1]
    edgelength.append(newedge+edgelength[i-1])
    i=i+1

#normalise:
edgelength=edgelength/edgelength[nedge-1]

#plot data:    
plt.plot(iedge,edgelength)

#plt.show()
saveplot = ('../../outputs/plots/cmd_snap'+str(snapnum))
plt.savefig(saveplot, bbox_inches='tight')
print "Graph saved at " + saveplot
plt.close()
