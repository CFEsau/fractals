#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import numpy as np
import sys
from sys import argv

snapnum, nmst, simnum, nran = argv[1:5] #snapnum passed from cdfplots.sh
#plot the cdf for this snapshot

snapnum=int(snapnum)
nmst=int(nmst)
nran=int(nran)

nedge=nmst-1

outpath = '../outputs/'
filepath = outpath+'runinv_k'+simnum+'/cluster_FoV5pc/'

#---------------------------
# CDF data for object stars
#---------------------------

objfile = 'MSTedges.dat'
objmacro = np.loadtxt(filepath+objfile)

#x-values:
objedgelength=[]
#Each edge length is recorded twice as the CDF here is a step function:
i=1
while i<nmst:
    objedgelength.append(objmacro[snapnum-1,i])
    objedgelength.append(objmacro[snapnum-1,i])
    i=i+1
    #snapnum is snapnum-1 because python is stupid and starts at 0.
    #Ditto 1st length being 'i' in array rather than 'i+1'.

#final x value is something arbitrary to get a horizontal line
#at the end of the CDF. Use 1.2*the longest edge length:
objedgelength.append(1.05*objedgelength[2*(nedge-1)])


#y-values:
probrange=np.linspace(0,1,num=nmst) #this gives nmst values from 0 to 1

#Again, CDF is a step function so need these values twice:
probability=[]
#1st y value is 0
probability.append(0)
#following y values are in the range 0>y>=1:
i=1
while i<nmst:
    probability.append(probrange[i])
    probability.append(probrange[i])
    i=i+1

#-----------------------
# CDFs for random stars
#-----------------------

j=1
while j<nran+1:
    
    filename = 'MSTedges_'+str(j)+'.dat'
    #data:
    macro = np.loadtxt(filepath+filename)

    
    #x-values:

    edgelength=[]

    #Each edge length is recorded twice as the CDF here is a step function:
    i=1
    while i<nmst:
        edgelength.append(macro[snapnum-1,i])
        edgelength.append(macro[snapnum-1,i])
        i=i+1
        #snapnum is snapnum-1 because python is stupid and starts at 0.
        #Ditto 1st length being 'i' in array rather than 'i+1'.

    #final x value is something arbitrary to get a horizontal line
    #at the end of the CDF. Use 1.2*the longest edge length:
    #if objedgelength[2*(nedge-1)]>edgelength[2*(nedge-1)]:
        #edgelength.append(1.05*objedgelength[2*(nedge-1)])
    #if objedgelength[2*(nedge-1)]>edgelength[2*(nedge-1)]:
    edgelength.append(1.05*edgelength[2*(nedge-1)])

    #y-values same as object stars calculated above

#-----------
# Plot CDFs
#-----------

    #my_dpi=96
    saveplot = ''
    plt.figure()
    plt.xlabel("Edge length (pc)")
    plt.ylabel("Probability")
    plt.title("Edge lengths of random stars")
    if objedgelength[2*(nedge-1)]>edgelength[2*(nedge-1)]:
        plt.xlim(xmax = 1.05*objedgelength[2*(nedge-1)], xmin = -0.01)
    if objedgelength[2*(nedge-1)]<edgelength[2*(nedge-1)]:
        plt.xlim(xmax = 1.05*edgelength[2*(nedge-1)], xmin = -0.01)
    plt.ylim(ymax = 1.005, ymin = 0)
    
    #plot data:
    #CDF for random stars
    plt.plot(edgelength,probability)
    #CDF for object stars
    plt.plot(objedgelength,probability,color='grey')#,linestyle='dashed')

    plt.annotate("snap "+str(snapnum), xy=(0.86,0.03), 
             xycoords='axes fraction', fontsize = 11)
    
    #plt.show()
    saveplot = (outpath+'plots/cdf_snap'+str(snapnum)+'_'+"%02d"%j+'.png')
    plt.savefig(saveplot, bbox_inches='tight')
    print "Graph saved at " + saveplot
    plt.close()
    
    j=j+1
