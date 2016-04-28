#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

print ("   Doing virial ratio...", sep="")

fbin, fdim_val, qvir_val = argv[1:4] #use when defining parameters inline

path = '../outputs'

duration = 10. #Duration of simulation (Myr)

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel("Qvir")
plt.title("Virial Ratio over Time")
plt.ylim(ymax = 1.0, ymin = 0.2)

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:
        ctype = 'all'
        
        filename = path + '/' + simname + '/energies.dat'
        macro = np.loadtxt(filename)
        nsnap =  macro[:,0]
        ekinetic = macro[:,1]
        epotential = macro[:,2]
        time = (nsnap/nsnap[-1])*duration
        qvir = ekinetic/-epotential
        
        #print ("   Doing ", simname, "...", sep="")
        plt.plot(time, qvir)
        
plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
             xycoords='axes fraction', fontsize = 10)
plt.annotate("fdim = " + fdim_val + "," + " qvir = " + qvir_val, 
             xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)

saveplot = path + '/plots/' + 'Qvir_comparek.pdf'
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print ("Saved in ", saveplot)
plt.close()
