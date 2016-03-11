#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

#print len(argv)
fbin, fdim_val, qvir_val = argv[1:4] #use when defining parameters inline

#path = str(raw_input('input the path to the simulation directory: '))
path = '../outputs'
#print path

#fdim_val = raw_input("fdim: ")
#qvir_val = raw_input("qvir: ")
#duration = float(input('input the duration of the simulation (Myr): '))
duration = 10

number = 0
while number <=0 or number > 5:
    print "Pick a graph:\n1 = energies\n2 = virial ratio\n3 = lagrangian radii\n4 = various lambda\n5 = all of the above"
    number = int(input('Pick a number: '))
    if  number <=0 or number > 5:
        print "\n Wrong number.\n"

for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        kval = simname.split("_")[1] #get k01, k02, etc
        filename1 = path + '/' + simname + '/macro'
        filename2 = path + '/' + simname + '/lambda'

        macro = np.loadtxt(filename1)
        nsnap =  macro[:,0]
        ekinetic = macro[:,1]
        epotential = macro[:,2]
        etotal = macro[:,3]
        lagrange = macro[:,4]
        time = (nsnap/nsnap[-1])*duration
        qvir = ekinetic/-epotential


        lambd = np.loadtxt(filename2)

        lam = lambd[:,1]
        lam_low = lambd[:,2]
        lam_up = lambd[:,3]
        yerr_low=lam-lam_low
        yerr_up=lam_up-lam
        lam_bar = lambd[:,4]
        lam_low_bar = lambd[:,5]
        lam_up_bar = lambd[:,6]
        yerr_low_bar = lam_bar-lam_low_bar
        yerr_up_bar = lam_up_bar-lam_bar
        lam_tilde = lambd[:,7]
        lam_low_tilde = lambd[:,8]
        lam_up_tilde = lambd[:,9]
        yerr_low_tilde = lam_tilde-lam_low_tilde
        yerr_up_tilde = lam_up_tilde-lam_tilde
        lam_star = lambd[:,10]
        lam_low_star = lambd[:,11]
        lam_up_star = lambd[:,12]
        yerr_low_star = lam_star-lam_low_star
        yerr_up_star = lam_up_star-lam_star
        gam = lambd[:,13]
        gam_low = lambd[:,14]
        gam_up = lambd[:,15]
        yerr_low_gam=gam-gam_low
        yerr_up_gam=gam_up-gam

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
            #plt.text(8.8,0.41e41,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,0.45e41,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_energies.pdf'
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
            plt.ylim(0,1)
            #plt.text(8.8,0.92,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,0.96,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_virial.pdf'
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
            #plt.text(8.8,3.7,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,3.86,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lagrangian.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

        elif number == 4:
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\Lambda$")
            plt.title("Mass segregation")
            #plt.plot(time, lam)
            plt.errorbar(time,lam,yerr=[yerr_low,yerr_up],errorevery=5)
            plt.ylim(0,15)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambda.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\overline{\Lambda}$")
            plt.title("Mass segregation")
            #plt.plot(time, lam_bar)
            plt.errorbar(time,lam_bar,yerr=[yerr_low_bar,yerr_up_bar],errorevery=5)
            plt.ylim(0,15)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambdabar.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\widetilde{\Lambda}$")
            plt.title("Mass segregation")
            #plt.plot(time, lam_tilde)
            plt.errorbar(time,lam_tilde,yerr=[yerr_low_tilde,yerr_up_tilde],errorevery=5)
            plt.ylim(0,25)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambdatilde.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\Lambda^\star$")
            plt.title("Mass segregation")
            #plt.plot(time, lam_star)
            plt.errorbar(time,lam_star,yerr=[yerr_low_star,yerr_up_star],errorevery=5)
            plt.ylim(0,15)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambdastar.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\Gamma$")
            plt.title("Mass segregation")
            #plt.plot(time, gam)
            plt.errorbar(time,gam,yerr=[yerr_low_gam,yerr_up_gam],errorevery=5)
            plt.ylim(0,25)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_gamma.pdf'
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
            #plt.text(8.8,0.41e41,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,0.45e41,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_energies.pdf'
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
            #plt.text(8.8,0.92,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,0.96,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_virial.pdf'
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
            #plt.text(8.8,3.7,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,3.86,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lagrangian.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
            
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Lambda")
            plt.title("Mass segregation")
            #plt.plot(time, lam)
            plt.errorbar(time,lam,yerr=[yerr_low,yerr_up],errorevery=5)
            plt.ylim(0,15)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambda.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\overline{\Lambda}$")
            plt.title("Mass segregation")
            #plt.plot(time, lam_bar)
            plt.errorbar(time,lam_bar,yerr=[yerr_low_bar,yerr_up_bar],errorevery=5)
            plt.ylim(0,15)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambdabar.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\widetilde{\Lambda}$")
            plt.title("Mass segregation")
            #plt.plot(time, lam_tilde)
            plt.errorbar(time,lam_tilde,yerr=[yerr_low_tilde,yerr_up_tilde],errorevery=5)
            plt.ylim(0,25)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambdatilde.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\Lambda^\star$")
            plt.title("Mass segregation")
            #plt.plot(time, lam_star)
            plt.errorbar(time,lam_star,yerr=[yerr_low_star,yerr_up_star],errorevery=5)
            plt.ylim(0,15)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_lambdastar.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()

            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\Gamma$")
            plt.title("Mass segregation")
            #plt.plot(time, gam)
            plt.errorbar(time,gam,yerr=[yerr_low_gam,yerr_up_gam],errorevery=5)
            plt.ylim(0,25)
            #plt.text(8.8,13.8,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,14.4,"qvir = " + str(qvir_val),fontsize=10)
            plt.annotate("fdim = " + str(fdim_val), xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + str(qvir_val), xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(str(kval), xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            saveplot = path + '/' + kval + '_gamma.pdf'
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print "Graph saved at " + saveplot
            plt.close()
