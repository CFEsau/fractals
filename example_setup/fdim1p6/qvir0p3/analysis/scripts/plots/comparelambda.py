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

#==============
#    Lambda
#==============

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\Lambda$")
plt.title("Mass segregation")
plt.ylim(0,15)

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        filename = path + '/' + simname + '/lambda'
        lambd = np.loadtxt(filename)
        nsnap =  lambd[:,0]
        time = (nsnap/nsnap[-1])*duration
        lam = lambd[:,1]

        print ("   Doing ", simname, "...", sep="")
        plt.plot(time, lam)

plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), xycoords='axes fraction', fontsize = 10)
plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
saveplot = path + "/comparelambda.pdf"
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print ("Saved in ", saveplot)
plt.close()
        

#==============
#  Lambda bar
#==============

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\overline{\Lambda}$")
plt.title("Mass segregation")
plt.ylim(0,15)

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        filename = path + '/' + simname + '/lambda'
        lambd = np.loadtxt(filename)
        nsnap =  lambd[:,0]
        time = (nsnap/nsnap[-1])*duration
        lam_bar = lambd[:,4]

        print ("   Doing ", simname, "...", sep="")
        plt.plot(time, lam_bar)

plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), xycoords='axes fraction', fontsize = 10)
plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
saveplot = path + "/comparelambdabar.pdf"
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print ("Saved in ", saveplot)
plt.close()
        

#==============
# Lambda tilde
#==============

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\widetilde{\Lambda}$")
plt.title("Mass segregation")
plt.ylim(0,25)

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        filename = path + '/' + simname + '/lambda'
        lambd = np.loadtxt(filename)
        nsnap =  lambd[:,0]
        time = (nsnap/nsnap[-1])*duration
        lam_tilde = lambd[:,7]

        print ("   Doing ", simname, "...", sep="")
        plt.plot(time, lam_tilde)

plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), xycoords='axes fraction', fontsize = 10)
plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
saveplot = path + "/comparelambdatilde.pdf"
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print ("Saved in ", saveplot)
plt.close()


#==============
# Lambda star
#==============

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\Lambda^\star$")
plt.title("Mass segregation")
plt.ylim(0,15)

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        filename = path + '/' + simname + '/lambda'
        lambd = np.loadtxt(filename)
        nsnap =  lambd[:,0]
        time = (nsnap/nsnap[-1])*duration
        lam_star = lambd[:,10]

        print ("   Doing ", simname, "...", sep="")
        plt.plot(time, lam_star)

plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), xycoords='axes fraction', fontsize = 10)
plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
saveplot = path + "/comparelambdastar.pdf"
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print ("Saved in ", saveplot)
plt.close()


#==============
#    Gamma
#==============

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\Gamma$")
plt.title("Mass segregation")
plt.ylim(0,15)

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        filename = path + '/' + simname + '/lambda'
        lambd = np.loadtxt(filename)
        nsnap =  lambd[:,0]
        time = (nsnap/nsnap[-1])*duration
        gam = lambd[:,13]

        print ("   Doing ", simname, "...", sep="")
        plt.plot(time, gam)

plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), xycoords='axes fraction', fontsize = 10)
plt.annotate('fdim = ' + fdim_val + ',' + ' qvir = ' + qvir_val, xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
saveplot = path + "/comparegamma.pdf"
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print ("Saved in ", saveplot)
plt.close()
