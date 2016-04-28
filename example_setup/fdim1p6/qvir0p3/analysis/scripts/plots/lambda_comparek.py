#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

fbin, fdim_val, qvir_val = argv[1:4] #use when defining parameters inline

path = '../outputs'

duration = 10. #Duration of simulation (Myr)

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''


#==============
#  Lambda bar
#==============

print ("   Doing lambda bar...", sep="")

for clustertype in os.listdir(path + '/runinv_k01/'):
    if 'cluster' in clustertype:

        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel(r"$\overline{\Lambda}$")
        plt.title("Mass segregation")
        plt.ylim(0,15)
        plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                     xycoords='axes fraction', fontsize = 10)
        plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, 
                     xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)

        for simname in os.listdir(path + '/'):
            if 'runinv' in simname:

                filename = (path + '/' + simname + '/' +
                            clustertype + '/lambda_xy')
                ctype = clustertype.split("_")[1] #get all, FoV, etc

                lambd = np.loadtxt(filename)
                nsnap =  lambd[:,0]
                time = (nsnap/nsnap[-1])*duration
                lam_bar = lambd[:,1]

                #print ("   Doing ", simname, "...", sep="")
                plt.plot(time, lam_bar)

        saveplot = path + '/plots/' + 'lambar_' + ctype + '_comparek.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print ("Saved in ", saveplot)
        plt.close()
        

#==============
# Lambda tilde
#==============
print ("")
print ("   Doing lambda tilde...", sep="")

for clustertype in os.listdir(path + '/runinv_k01/'):
    if 'cluster' in clustertype:

        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel(r"$\widetilde{\Lambda}$")
        plt.title("Mass segregation")
        plt.ylim(0,25)
        plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                     xycoords='axes fraction', fontsize = 10)
        plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, 
                     xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)

        for simname in os.listdir(path + '/'):
            if 'runinv' in simname:

                filename = (path + '/' + simname + '/' +
                            clustertype + '/lambda_xy')
                ctype = clustertype.split("_")[1] #get all, FoV, etc

                lambd = np.loadtxt(filename)
                nsnap =  lambd[:,0]
                time = (nsnap/nsnap[-1])*duration
                lam_til = lambd[:,4]

                #print ("   Doing ", simname, "...", sep="")
                plt.plot(time, lam_til)

        saveplot = path + '/plots/' + 'lamtil_' + ctype + '_comparek.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print ("Saved in ", saveplot)
        plt.close()


#==============
# Lambda star
#==============
print ("")
print ("   Doing lambda star...", sep="")

for clustertype in os.listdir(path + '/runinv_k01/'):
    if 'cluster' in clustertype:

        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel(r"$\Lambda^\star$")
        plt.title("Mass segregation")
        plt.ylim(0,15)
        plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                     xycoords='axes fraction', fontsize = 10)
        plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, 
                     xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)

        for simname in os.listdir(path + '/'):
            if 'runinv' in simname:

                filename = (path + '/' + simname + '/' +
                            clustertype + '/lambda_xy')
                ctype = clustertype.split("_")[1] #get all, FoV, etc

                lambd = np.loadtxt(filename)
                nsnap =  lambd[:,0]
                time = (nsnap/nsnap[-1])*duration
                lam_star = lambd[:,7]

                #print ("   Doing ", simname, "...", sep="")
                plt.plot(time, lam_star)

        saveplot = path + '/plots/' + 'lamstar_' + ctype + '_comparek.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print ("Saved in ", saveplot)
        plt.close()


#==============
#    Gamma
#==============
print ("")
print ("   Doing gamma...", sep="")

for clustertype in os.listdir(path + '/runinv_k01/'):
    if 'cluster' in clustertype:

        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel(r"$\Gamma$")
        plt.title("Mass segregation")
        plt.ylim(0,25)
        plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                     xycoords='axes fraction', fontsize = 10)
        plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, 
                     xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)

        for simname in os.listdir(path + '/'):
            if 'runinv' in simname:

                filename = (path + '/' + simname + '/' +
                            clustertype + '/lambda_xy')
                ctype = clustertype.split("_")[1] #get all, FoV, etc

                lambd = np.loadtxt(filename)
                nsnap =  lambd[:,0]
                time = (nsnap/nsnap[-1])*duration
                gam = lambd[:,10]

                #print ("   Doing ", simname, "...", sep="")
                plt.plot(time, gam)

        saveplot = path + '/plots/' + 'gamma_' + ctype + '_comparek.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print ("Saved in ", saveplot)
        plt.close()
