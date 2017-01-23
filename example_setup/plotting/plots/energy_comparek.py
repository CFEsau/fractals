#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

print ("   Doing energies...", sep="")

fbin, fdim_val, qvir_val, outpath = argv[1:5] #use when defining parameters inline

path = outpath + '/outputs/'

duration = 10. #Duration of simulation (Myr)

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel("Energy (J)")
plt.title("Energies over Time")
plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
#plt.text(8.4,0.45e41,"qvir =" + str(qvir_val))
#plt.text(8.4,0.4e41,"fdim =" + str(fdim_val))

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:
        ctype = 'all'

        filename = path + '/' + simname + '/energies.dat'
        macro = np.loadtxt(filename)
        nsnap =  macro[:,0]
        ekinetic = macro[:,1]
        epotential = macro[:,2]
        etotal = macro[:,3]
        time = (nsnap/nsnap[-1])*duration

        #print ("   Doing ", simname, "...", sep="")
        plt.plot(time, ekinetic, label = "Kinetic Energy")
        plt.plot(time, epotential, label = "Potential Energy")
        plt.plot(time, etotal, label = "Total Energy")

#plt.legend(loc=2)
plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), xycoords='axes fraction', fontsize = 10)
plt.annotate("fdim = " + fdim_val + "," + " qvir = " + qvir_val, xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)

saveplot = path + '/plots/' + 'E_allk.pdf'
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print ("Saved in ", saveplot)
plt.close()
