#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

fbin, fdim_val, qvir_val = argv[1:4]

path = '../outputs'

duration = 10 #Duration of simulation (Myr)


for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        kval = simname.split("_")[1] #get k01, k02, etc
        filename = path + '/' + simname + '/energies.dat'
        ctype = 'all'

        macro = np.loadtxt(filename)
        nsnap =  macro[:,0]
        ekinetic = macro[:,1]
        epotential = macro[:,2]
        time = (nsnap/nsnap[-1])*duration
        qvir = ekinetic/-epotential


        my_dpi=96
        #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)

        saveplot = ''
        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel("Qvir")
        plt.title("Virial Ratio over Time")
        plt.plot(time, qvir)
        plt.ylim(0,1)
        #plt.text(8.8,0.92,"fdim = " + str(fdim_val),fontsize=10)
        #plt.text(8.8,0.96,"qvir = " + str(qvir_val),fontsize=10)
        plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                     xycoords='axes fraction', horizontalalignment='right',
                     verticalalignment='bottom', fontsize=10)
        plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                     xycoords='axes fraction', horizontalalignment='right',
                     verticalalignment='bottom', fontsize=10)
        plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                     xycoords='axes fraction', horizontalalignment='right',
                     verticalalignment='bottom', fontsize=10)
        plt.annotate(str(kval), xy=(0.99, -0.08),
                     xycoords='axes fraction', horizontalalignment='right',
                     verticalalignment='bottom', fontsize=10)

        saveplot = path + '/plots/' + kval + '_' + ctype + '_Qvir.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at " + saveplot
        plt.close()
