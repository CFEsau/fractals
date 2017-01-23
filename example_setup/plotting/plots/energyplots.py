#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

fbin, fdim_val, qvir_val, outpath = argv[1:5]

path = outpath + '/outputs/'

duration = 10. #Duration of simulation (Myr)


for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        kval = simname.split("_")[1] #get k01, k02, etc
        ctype = 'all'
        filename = path + simname + '/energies.dat'

        macro = np.loadtxt(filename)
        nsnap =  macro[:,0]
        ekinetic = macro[:,1]
        epotential = macro[:,2]
        etotal = macro[:,3]
        time = (nsnap/nsnap[-1])*duration

        my_dpi=96
        #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
        
        saveplot = ''
        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel("Energy (J)")
        plt.title("Energies over Time")
        plt.plot(time, ekinetic, label = 'Kinetic Energy')
        plt.plot(time, epotential, label = 'Potential Energy')
        plt.plot(time, etotal, label = 'Total Energy')
        plt.legend(loc=2,fontsize=11)
        plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
        #plt.text(8.8,0.41e41,"fdim = " + str(fdim_val),fontsize=10)
        #plt.text(8.8,0.45e41,"qvir = " + str(qvir_val),fontsize=10)
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

        saveplot = path + 'plots/' + kval + '_' + ctype + '_E.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at " + saveplot
        plt.close()
