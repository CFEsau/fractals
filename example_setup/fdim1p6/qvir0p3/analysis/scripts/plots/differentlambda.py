#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

fbin, fdim_val, qvir_val = argv[1:4] #use when defining parameters inline

#fdim_val = str(input('Enter fdim: '))
#qvir_val = str(input('Enter qvir: '))
#qparam1,qparam2=qvir_val.split(".")
#qparam_str=qparam1 + 'p' + qparam2
#path = 'qvir' + qparam_str + '/outputs/'

#path = str(raw_input('input the path to the simulation directory: '))
path = '../outputs'
#print (path)

#duration = float(input('Input the duration of the simulation (Myr): '))
duration = 10.

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

#===================================================
#  Comparing different measures of mass segregation
#===================================================


for simname in os.listdir(path + '/'):
    if 'runinv' in simname:
        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel(r"$\Lambda$")
        plt.title("Mass segregation")
        plt.ylim(0,25)

        kval = simname.split("_")[1] #get k01, k02, etc
        filename = path + '/' + simname + '/lambda'
        lambd = np.loadtxt(filename)
        nsnap =  lambd[:,0]
        time = (nsnap/nsnap[-1])*duration
        lam_bar = lambd[:,4]
        lam_tilde = lambd[:,7]
        lam_star = lambd[:,10]
        gam = lambd[:,13]

        print ("Doing ", simname, "...", sep="")
        plt.plot(time, lam_bar, label = r"$\overline{\Lambda}$")
        plt.plot(time, lam_tilde, label = r"$\widetilde{\Lambda}$")
        plt.plot(time, lam_star, label = r"$\Lambda^\star$")
        plt.plot(time, gam, label = r"$\Gamma$")


        plt.legend(fontsize=11)
        plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), xycoords='axes fraction', fontsize = 10)
        plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
        saveplot = path + '/' + kval + "_differentlambda.pdf"
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print ("   Graph saved at ", saveplot)
        plt.close()
