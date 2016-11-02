#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

print ("   Doing half-mass radius...", sep="")

fbin, fdim_val, qvir_val = argv[1:4] #use when defining parameters inline

path = '../outputs'

duration = 10. #Duration of simulation (Myr)

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

for clustype in os.listdir(path + '/runinv_k01/'):
    if 'cluster_' in clustype:

        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel("Lagrangian Radius (pc)")
        plt.title("Lagrangian Radii over Time")
        plt.ylim(ymax = 7, ymin = 0)
        label = '0.5'
        #plt.legend(loc=2)
        plt.annotate("r = " + label,xy=(0.1,0.9), 
                     xycoords='axes fraction', fontsize = 12)
        plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                     xycoords='axes fraction', fontsize = 10)
        plt.annotate("fdim = " + fdim_val + "," + " qvir = " + qvir_val, 
                     xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)

        for simname in os.listdir(path + '/'):
            if 'runinv' in simname:

                ctype = clustype.split("_")[1] #get all, FoV, etc
                filename = (path + '/' + simname + '/' + clustype + 
                            '/c_of_m_3D.dat')

                macro = np.loadtxt(filename)
                nsnap =  macro[:,0]
                lagrange = macro[:,4]
                time = (nsnap/nsnap[-1])*duration
                #print ("   Doing ", simname, "...", sep="")
                plt.plot(time, lagrange)

        saveplot = path + '/plots/' + 'Rh_' + ctype + '_comparek.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print ("Saved in ", saveplot)
        plt.close()
