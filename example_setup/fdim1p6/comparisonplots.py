#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

nsim = int(input('Enter the number of simulations: '))
path1 = str('qvir0p3/kdum')
path2 = str('/analysis/outputs/runinv_0100')

duration = float(input('Input the duration of the simulation (Myr): '))

number = 0
while number <=0 or number > 4:
        print "Pick a graph:\n1 = energies\n2 = virial ratio\n3 = lagrangian radii\n4 = all of the above"
        number = int(input('Pick a number: '))
        if  number <=0 or number > 4:
                print "\n Wrong number.\n"

my_dpi=96
plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

if number == 1:
        plt.xlabel('Time (Myr)')
        plt.ylabel('Energy (J)')
        plt.title('Energies over Time')
        plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
        for k in range(1,nsim+1):
                filename = path1 + str(k) + path2 + "/macro"
                macro = np.loadtxt(filename)
                nsnap =  macro[:,0]
                ekinetic = macro[:,1]
                epotential = macro[:,2]
                etotal = macro[:,3]
                lagrange = macro[:,4]
                time = (nsnap/nsnap[-1])*duration
                qvir = ekinetic/-epotential
                print "  Doing",filename
                plt.plot(time, ekinetic, label = 'Kinetic Energy')
                plt.plot(time, epotential, label = 'Potential Energy')
                plt.plot(time, etotal, label = 'Total Energy')
        #plt.legend(loc=2)
        saveplot = 'compareenergies.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at ", saveplot
        
elif number == 2:
        plt.xlabel('Time (Myr)')
        plt.ylabel('Qvir')
        plt.title('Virial Ratio over Time')
        plt.ylim(ymax = 1.0, ymin = 0.2)
        for k in range(1,nsim+1):
                filename = path1 + str(k) + path2 + "/macro"
                macro = np.loadtxt(filename)
                nsnap =  macro[:,0]
                ekinetic = macro[:,1]
                epotential = macro[:,2]
                etotal = macro[:,3]
                lagrange = macro[:,4]
                time = (nsnap/nsnap[-1])*duration
                qvir = ekinetic/-epotential
                print "  Doing",filename       
                plt.plot(time, qvir)
        saveplot = 'comparevirial.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at ", saveplot

elif number == 3:
        plt.xlabel('Time (Myr)')
        plt.ylabel('Lagrangian Radius (pc)')
        plt.title('Lagrangian Radii over Time')
        plt.ylim(ymax = 7, ymin = 0)
        label = '0.5'
        for k in range(1,nsim+1):
                filename = path1 + str(k) + path2 + "/macro"
                macro = np.loadtxt(filename)
                nsnap =  macro[:,0]
                ekinetic = macro[:,1]
                epotential = macro[:,2]
                etotal = macro[:,3]
                lagrange = macro[:,4]
                time = (nsnap/nsnap[-1])*duration
                qvir = ekinetic/-epotential
                print "  Doing",filename
                plt.plot(time, lagrange, label = 'r = ' + label)
        #plt.legend(loc=2)
        plt.annotate('r = ' + label,xy=(0.1,0.9),xycoords='axes fraction', fontsize = 12)
        saveplot = 'comparelagrangian.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at ", saveplot
        
elif number == 4:
        plt.figure(1)
        plt.xlabel('Time (Myr)')
        plt.ylabel('Energy (J)')
        plt.title('Energies over Time')
        plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
        for k in range(1,nsim+1):
                filename = path1 + str(k) + path2 + "/macro"
                macro = np.loadtxt(filename)
                nsnap =  macro[:,0]
                ekinetic = macro[:,1]
                epotential = macro[:,2]
                etotal = macro[:,3]
                lagrange = macro[:,4]
                time = (nsnap/nsnap[-1])*duration
                qvir = ekinetic/-epotential
                print "  Doing",filename
                plt.plot(time, ekinetic, label = 'Kinetic Energy')
                plt.plot(time, epotential, label = 'Potential Energy')
                plt.plot(time, etotal, label = 'Total Energy')
        #plt.legend(loc=2)
        saveplot = 'compareenergies.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at ", saveplot
        
        plt.figure(2)
        plt.xlabel('Time (Myr)')
        plt.ylabel('Qvir')
        plt.title('Virial Ratio over Time')
        plt.ylim(ymax=1.0, ymin = 0.2)
        for k in range(1,nsim+1):
                filename = path1 + str(k) + path2 + "/macro"
                macro = np.loadtxt(filename)
                nsnap =  macro[:,0]
                ekinetic = macro[:,1]
                epotential = macro[:,2]
                etotal = macro[:,3]
                lagrange = macro[:,4]
                time = (nsnap/nsnap[-1])*duration
                qvir = ekinetic/-epotential
                print "  Doing",filename       
                plt.plot(time, qvir)
        saveplot = 'comparevirial.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at ", saveplot
        
        plt.figure(3)
        plt.xlabel('Time (Myr)')
        plt.ylabel('Lagrangian Radius (pc)')
        plt.title('Lagrangian Radii over Time')
        plt.ylim(ymax = 7, ymin = 0)
        label = '0.5'
        for k in range(1,nsim+1):
                filename = path1 + str(k) + path2 + "/macro"
                macro = np.loadtxt(filename)
                nsnap =  macro[:,0]
                ekinetic = macro[:,1]
                epotential = macro[:,2]
                etotal = macro[:,3]
                lagrange = macro[:,4]
                time = (nsnap/nsnap[-1])*duration
                qvir = ekinetic/-epotential
                print "  Doing",filename
                plt.plot(time, lagrange, label = 'r = ' + label)
        #plt.legend(loc=2)
        plt.annotate('r = ' + label,xy=(0.1,0.9),xycoords='axes fraction', fontsize = 12)
        saveplot = 'comparelagrangian.pdf'
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        print "Graph saved at ", saveplot
        
print "     ...done"
#plt.tight_layout()
#plt.savefig(path + saveplot, bbox_inches='tight')

#print "Graph saved at " + path + saveplot
