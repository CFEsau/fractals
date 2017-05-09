#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import plotconfig
    
def printenergy():
    print ("")
    print ("------------")
    print ("-- Energy --")
    print ("------------")
    
def plotenergy():
    duration = 10. #Duration of simulation (Myr)
    for simname in os.listdir(plotconfig.outpath + '/'):
        if 'runinv' in simname:
            
            kval = simname.split("_")[1] #get k01, k02, etc
            ctype = 'all'
            
            filename = plotconfig.outpath + '/' + simname + '/energies.dat'
            
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            etotal = macro[:,3]
            time = (nsnap/nsnap[-1])*duration
            
            my_dpi=96
            #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 is 8.3 x 11.7 (portrait)
            
            saveplot = ''
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
            plt.annotate("fdim = " + plotconfig.fdim, xy=(0.99, 0.96),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("qvir = " + plotconfig.qvir, xy=(0.99, 0.92),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + plotconfig.fbin, xy=(0., -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(kval, xy=(0.99, -0.08),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            
            saveplot = (plotconfig.outpath + '/plots/' +
                        kval + '_' + ctype + '_E.pdf')
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print ("    ",kval,"_",ctype,"_E.pdf",sep="")
            plt.close()

def energy_k():
    print ("   Comparing across k...")
    
    duration = 10. #Duration of simulation (Myr)
    my_dpi=96
    #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 landscape
    saveplot = ''
    
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Energy (J)")
    plt.title("Energies over Time")
    plt.ylim(ymax = 0.5e41, ymin = -0.6e41)
    #plt.text(8.4,0.45e41,"qvir =" + str(qvir_val))
    #plt.text(8.4,0.4e41,"fdim =" + str(fdim_val))
    
    for simname in os.listdir(plotconfig.outpath + '/'):
        if 'runinv' in simname:
            ctype = 'all'
            filename = plotconfig.outpath + '/' + simname + '/energies.dat'
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
    plt.annotate("fbin = " + plotconfig.fbin, xy=(0.95,-0.09),
                 xycoords='axes fraction', fontsize = 10)
    plt.annotate("fdim = "+plotconfig.fdim+", qvir = "+plotconfig.qvir,
                 xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
    
    saveplot = plotconfig.outpath + '/plots/E_comparek.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("         E_comparek.pdf")
    print ("")
    plt.close()
