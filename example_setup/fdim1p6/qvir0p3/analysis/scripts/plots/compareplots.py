#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

fbin, fdim_val, qvir_val = argv[1:4] #use when defining parameters inline

#fdim = str(input('Enter fdim: '))
#qparam = str(input('Enter qvir: '))
qparam1,qparam2=qparam.split(".")
qparam_str=qparam1 + 'p' + qparam2
path = 'qvir' + qparam_str + '/outputs/'

#duration = float(input('Input the duration of the simulation (Myr): '))
duration = 10.

number = 0
while number <=0 or number > 5:
        print ("Pick a graph:\n1 = energies\n2 = virial ratio\n3 = lagrangian radii\n4 = lambda\n5 = all of the above")
        number = int(input("Pick a number: "))
        if  number <=0 or number > 5:
                print ("\n Wrong number.\n")

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

if number == 1:
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy (J)")
    plt.title("Energies over Time")
    plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
    for kdir in os.listdir(path):

        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            etotal = macro[:,3]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, ekinetic, label = 'Kinetic Energy')
            plt.plot(time, epotential, label = 'Potential Energy')
            plt.plot(time, etotal, label = 'Total Energy')

    #plt.legend(loc=2)
    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'compareenergies_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()

elif number == 2:
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Qvir")
    plt.title("Virial Ratio over Time")
    plt.ylim(ymax = 1.0, ymin = 0.2)
    for kdir in os.listdir(path):

        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            etotal = macro[:,3]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, qvir)

    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'comparevirial_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()

elif number == 3:
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Lagrangian Radius (pc)")
    plt.title("Lagrangian Radii over Time")
    plt.ylim(ymax = 7, ymin = 0)
    label = '0.5'
    for kdir in os.listdir(path):

        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            etotal = macro[:,3]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, lagrange, label = 'r = ' + label)

    #plt.legend(loc=2)
    plt.annotate('r = ' + label,xy=(0.1,0.9),xycoords='axes fraction', fontsize = 12)
    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'comparelagrangian_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()

elif number == 4:
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Lambda")
    plt.title("Mass segregation")
    plt.ylim(0,15)
    for kdir in os.listdir(path):

        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            time = (nsnap/nsnap[-1])*duration
            lam = macro[:,5]
            lam_low = macro[:,6]
            lam_up = macro[:,7]
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, lam)

    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'comparelambda_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()

elif number == 5:
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy (J)")
    plt.title("Energies over Time")
    plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
    for kdir in os.listdir(path):

        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            etotal = macro[:,3]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, ekinetic, label = 'Kinetic Energy')
            plt.plot(time, epotential, label = 'Potential Energy')
            plt.plot(time, etotal, label = 'Total Energy')

    #plt.legend(loc=2)
    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'compareenergies_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()
        
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Qvir")
    plt.title("Virial Ratio over Time")
    plt.ylim(ymax=1.0, ymin = 0.2)
    for kdir in os.listdir(path):

        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            etotal = macro[:,3]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, qvir)

    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'comparevirial_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()

    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Lagrangian Radius (pc)")
    plt.title("Lagrangian Radii over Time")
    plt.ylim(ymax = 7, ymin = 0)
    label = '0.5'
    for kdir in os.listdir(path):

        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            etotal = macro[:,3]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, lagrange, label = 'r = ' + label)

    #plt.legend(loc=2)
    plt.annotate('r = ' + label,xy=(0.1,0.9),xycoords='axes fraction', fontsize = 12)
    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'comparelagrangian_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()

    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Lambda")
    plt.title("Mass segregation")
    plt.ylim(0,15)
    for kdir in os.listdir(path):
        if 'runinv' in kdir:
            filename = path + kdir + '/macro'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            time = (nsnap/nsnap[-1])*duration
            lam = macro[:,5]
            lam_low = macro[:,6]
            lam_up = macro[:,7]
            print ("   Doing ", kdir, "...", sep="")
            plt.plot(time, lam)
    plt.annotate('fdim = ' + fdim + ',' + ' qvir = ' + qparam, xy=(-0.1,-0.1), xycoords='axes fraction', fontsize = 12)
    saveplot = 'comparelambda_qvir' + qparam_str + '.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("Graph saved at ", saveplot)
    plt.close()
        
print ("\n                ...done\n")
