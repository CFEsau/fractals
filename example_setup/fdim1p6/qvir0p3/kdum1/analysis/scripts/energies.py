#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

path = str(raw_input('input the path to the simulation directory: '))
print path
filename = path + "/macro"

duration = float(input('input the duration of the simulation (Myr): '))

macro = np.loadtxt(filename)
nsnap =  macro[:,0]
ekinetic = macro[:,1]
epotential = macro[:,2]
etotal = macro[:,3]
lagrange = macro[:,4]
time = (nsnap/nsnap[-1])*duration
qvir = ekinetic/-epotential

number = 0
while number <=0 or number > 4:
    print "Pick a graph:\n1 = energies\n2 = virial ratio\n3 = lagrangian radii\n4 = all of the above"
    number = int(input('Pick a number: '))
    if  number <=0 or number > 4:
        print "\n Wrong number.\n"
    
my_dpi=96
plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)

plt.xlabel('Time (Myr)')
saveplot = ''
if number == 1:
    plt.ylabel('Energy (J)')
    plt.title('Energies over Time')
    plt.plot(time, ekinetic, label = 'Kinetic Energy')
    plt.plot(time, epotential, label = 'Potential Energy')
    plt.plot(time, etotal, label = 'Total Energy')
    plt.legend(loc=2)
    saveplot = '/energies.pdf'
    plt.tight_layout()
    plt.savefig(path + saveplot, bbox_inches='tight')
    print "Graph saved at " + path + saveplot
elif number == 2:
    plt.ylabel('Qvir')
    plt.title('Virial Ratio over Time')
    plt.plot(time, qvir)
    saveplot = '/virial.pdf'
    plt.tight_layout()
    plt.savefig(path + saveplot, bbox_inches='tight')
    print "Graph saved at " + path + saveplot
elif number == 3:
    plt.ylabel('Lagrangian Radius (pc)')
    plt.title('Lagrangian Radii over Time')
    label = '0.5'
    plt.plot(time, lagrange, label = 'r = ' + label)
    plt.legend(loc=2)
    saveplot = '/lagrangian.pdf'
    plt.tight_layout()
    plt.savefig(path + saveplot, bbox_inches='tight')
    print "Graph saved at " + path + saveplot
elif number == 4:
    plt.figure(1)
    plt.ylabel('Energy (J)')
    plt.title('Energies over Time')
    plt.plot(time, ekinetic, label = 'Kinetic Energy')
    plt.plot(time, epotential, label = 'Potential Energy')
    plt.plot(time, etotal, label = 'Total Energy')
    plt.legend(loc=2)
    saveplot = '/energies.pdf'
    plt.tight_layout()
    plt.savefig(path + saveplot, bbox_inches='tight')
    print "Graph saved at " + path + saveplot

    plt.figure(2)
    plt.ylabel('Qvir')
    plt.title('Virial Ratio over Time')
    plt.plot(time, qvir)
    saveplot = '/virial.pdf'
    plt.tight_layout()
    plt.savefig(path + saveplot, bbox_inches='tight')
    print "Graph saved at " + path + saveplot

    plt.figure(3)
    plt.ylabel('Lagrangian Radius (pc)')
    plt.title('Lagrangian Radii over Time')
    label = '0.5'
    plt.plot(time, lagrange, label = 'r = ' + label)
    plt.legend(loc=2)
    saveplot = '/lagrangian.pdf'
    plt.tight_layout()
    plt.savefig(path + saveplot, bbox_inches='tight')
    print "Graph saved at " + path + saveplot

#plt.tight_layout()
#plt.savefig(path + saveplot, bbox_inches='tight')

#print "Graph saved at " + path + saveplot
