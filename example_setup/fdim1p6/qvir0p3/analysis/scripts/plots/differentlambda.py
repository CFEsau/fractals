#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

fbin, fdim_val, qvir_val = argv[1:4] #use when defining parameters inline

path = '../outputs'

duration = 10. #duration of simulation (Myr)

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

#===================================================
#  Comparing different measures of mass segregation
#===================================================


for simname in os.listdir(path + '/'):
    if 'runinv' in simname:
        for clustype in os.listdir(path + '/' + simname + '/'):
            if 'cluster_' in clustype:

                kval = simname.split("_")[1] #get k01, k02, etc
                filename = path + '/' + simname + '/' + clustype + '/lambda_xy'
                lambd = np.loadtxt(filename)
                nsnap =  lambd[:,0]
                time = (nsnap/nsnap[-1])*duration
                lam_bar = lambd[:,1]
                lam_tilde = lambd[:,4]
                lam_star = lambd[:,7]
                gam = lambd[:,10]

                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(r"$\Lambda$")
                plt.title("Mass segregation")
                plt.ylim(0,25)

                plt.plot(time, lam_bar, label = r"$\overline{\Lambda}$")
                plt.plot(time, lam_tilde, label = r"$\widetilde{\Lambda}$")
                plt.plot(time, lam_star, label = r"$\Lambda^\star$")
                plt.plot(time, gam, label = r"$\Gamma$")


                plt.legend(fontsize=10)
                plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                             xycoords='axes fraction', fontsize = 10)
                plt.annotate('fdim = ' + fdim_val + ', qvir = ' + qvir_val, 
                             xy=(-0.06,-0.09), xycoords='axes fraction', 
                             fontsize = 10)

                clusterstring = clustype.split("_")[1]

                saveplot = (path + '/plots/' + 
                            kval + '_' + clusterstring + '_lam.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                print ("   Graph saved at ", saveplot)
                plt.close()
