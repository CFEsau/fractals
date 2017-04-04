#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import plotconfig
    
def printvirial():
    print ("")
    print ("------------------")
    print ("-- Virial ratio --")
    print ("------------------")

def plotvirial():
    duration = 10. #Duration of simulation (Myr)
    for simname in os.listdir(plotconfig.outpath + '/'):
        if 'runinv' in simname:
            
            kval = simname.split("_")[1] #get k01, k02, etc
            filename = plotconfig.outpath + '/' + simname + '/energies.dat'
            ctype = 'all'
            
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            
            my_dpi=96
            #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 landscape
            
            saveplot = ''
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel("Qvir")
            plt.title("Virial Ratio over Time")
            plt.plot(time, qvir)
            plt.ylim(0,1)
            #plt.text(8.8,0.92,"fdim = " + str(fdim_val),fontsize=10)
            #plt.text(8.8,0.96,"qvir = " + str(qvir_val),fontsize=10)
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
                        kval + '_' + ctype + '_Qvir.pdf')
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print ("    ",kval,"_",ctype,"_Qvir.pdf",sep="")
            plt.close()
            
def virial_k():
    print ("   Doing virial ratio...")

    duration = 10. #Duration of simulation (Myr)
    my_dpi=96
    #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 is 8.3 x 11.7 (portrait)
    saveplot = ''
    
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Qvir")
    plt.title("Virial Ratio over Time")
    plt.ylim(ymax = 1.0, ymin = 0.2)
    
    for simname in os.listdir(plotconfig.outpath + '/'):
        if 'runinv' in simname:
            ctype = 'all' 
            filename = plotconfig.outpath + '/' + simname + '/energies.dat'
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            ekinetic = macro[:,1]
            epotential = macro[:,2]
            time = (nsnap/nsnap[-1])*duration
            qvir = ekinetic/-epotential
            
            #print ("   Doing ", simname, "...", sep="")
            plt.plot(time, qvir)
            
    plt.annotate("fbin = " + plotconfig.fbin, xy=(0.95,-0.09), 
                 xycoords='axes fraction', fontsize = 10)
    plt.annotate("fdim = "+plotconfig.fdim+", qvir = "+plotconfig.qvir,
                 xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
    
    saveplot = plotconfig.outpath + '/plots/Qvir_comparek.pdf'
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ("         Qvir_comparek.pdf")
    print ("")
    plt.close()
