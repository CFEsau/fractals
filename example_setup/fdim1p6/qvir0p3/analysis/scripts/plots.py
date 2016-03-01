#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

path = str(raw_input('input the path to the simulation directory: '))
print path

fdim_val = raw_input("fdim: ")
qvir_val = raw_input("qvir: ")
#duration = float(input('input the duration of the simulation (Myr): '))
duration = 10

number = 0
while number <=0 or number > 5:
    print "Pick a graph:\n1 = energies\n2 = virial ratio\n3 = lagrangian radii\n4 = lambda\n5 = all of the above"
    number = int(input('Pick a number: '))
    if  number <=0 or number > 5:
        print "\n Wrong number.\n"

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        filename = path + '/' + simname + '/macro'

        macro = np.loadtxt(filename)
        nsnap =  macro[:,0]
        ekinetic = macro[:,1]
        epotential = macro[:,2]
        etotal = macro[:,3]
        lagrange = macro[:,4]
        time = (nsnap/nsnap[-1])*duration
        qvir = ekinetic/-epotential
        lam = macro[:,5]
        lam_low = macro[:,6]
        lam_up = macro[:,7]
        
        my_dpi=96
        #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
        
        saveplot = ''
        if number == 1:
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Energy (J)")
            plt.title("Energies over Time")
            plt.plot(time, ekinetic, label = 'Kinetic Energy')
            plt.plot(time, epotential, label = 'Potential Energy')
            plt.plot(time, etotal, label = 'Total Energy')
            plt.legend(loc=2,fontsize=11)
            plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
            plt.text(8.8,0.45e41,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,0.4e41,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/energies.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
        elif number == 2:
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Qvir")
            plt.title("Virial Ratio over Time")
            plt.plot(time, qvir)
            plt.ylim((0,1))
            plt.text(8.8,0.96,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,0.92,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/virial.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
        elif number == 3:
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Lagrangian Radius (pc)")
            plt.title("Lagrangian Radii over Time")
            label = '0.5'
            plt.plot(time, lagrange, label = 'r = ' + label)
            plt.legend(loc=2,fontsize=11)
            plt.ylim(0,4)
            plt.text(8.8,3.86,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,3.7,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/lagrangian.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
        elif number == 4:
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Lambda")
            plt.title("Mass segregation")
            plt.plot(time, lam)
            plt.ylim(0,15)
            plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/lambda.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
            
        elif number == 5:
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Energy (J)")
            plt.title("Energies over Time")
            plt.plot(time, ekinetic, label = 'Kinetic Energy')
            plt.plot(time, epotential, label = 'Potential Energy')
            plt.plot(time, etotal, label = 'Total Energy')
            plt.legend(loc=2,fontsize=11)
            plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
            plt.text(8.8,0.45e41,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,0.4e41,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/energies.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
            
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Qvir")
            plt.title("Virial Ratio over Time")
            plt.plot(time, qvir)
            plt.ylim(0,1)
            plt.text(8.8,0.96,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,0.92,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/virial.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
            
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Lagrangian Radius (pc)")
            plt.title("Lagrangian Radii over Time")
            label = '0.5'
            plt.plot(time, lagrange, label = 'r = ' + label)
            plt.legend(loc=2,fontsize=11)
            plt.ylim(0,4)
            plt.text(8.8,3.86,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,3.7,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/lagrangian.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
            
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Lambda")
            plt.title("Mass segregation")
            plt.plot(time, lam)
            plt.ylim((0,15))
            plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            saveplot = path + '/' + simname + '/lambda.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
            
#plt.tight_layout()
#plt.savefig(saveplot, bbox_inches='tight')
            
#print "Graph saved at " + saveplot
