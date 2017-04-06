#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import plotconfig
    
def printrhalfm():
    print ("")
    print ("----------------------")
    print ("-- Half-mass radius --")
    print ("----------------------")

def plotrhalfm(thiscluster):
    duration = 10. #Duration of simulation (Myr)
    
    for simname in os.listdir(plotconfig.outpath + '/'):
        if 'runinv' in simname:
            kval = simname.split("_")[1] #get k01, k02, etc
            
            #see whether cluster type has been done for this k:
            #if os.path.exists(plotconfig.outpath+'/'+simname+'/'+thiscluster):
            ctype = thiscluster.split("_")[1] #get all, FoV, etc
            filename = (plotconfig.outpath+'/'+simname+'/' + 
                        thiscluster + '/c_of_m_3D.dat')
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            
            my_dpi=96
            #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 landscape
            
            saveplot = ''
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
                        kval + '_' + ctype + '_Rh.pdf')
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            print ("    ",kval,"_",ctype,"_Rh.pdf",sep="")
            plt.close()
                
def rhalfm_k(thiscluster):
    print ("   Comparing across k...")
    my_dpi=96
    duration = 10. #Duration of simulation (Myr)
    plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 landscape
    saveplot = ''
    
    plt.figure()
    plt.xlabel("Time (Myr)")
    plt.ylabel("Lagrangian Radius (pc)")
    plt.title("Lagrangian Radii over Time")
    plt.ylim(ymax = 7, ymin = 0)
    label = '0.5'
    
    #plt.legend(loc=2)
    plt.annotate("r = " + label,xy=(0.1,0.9), 
                 xycoords='axes fraction', fontsize = 12)
    plt.annotate("fbin = " + plotconfig.fbin, xy=(0.95,-0.09), 
                 xycoords='axes fraction', fontsize = 10)
    plt.annotate("fdim = " + plotconfig.fdim + ", qvir = " + plotconfig.qvir,
                 xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
    
    for simname in os.listdir(plotconfig.outpath + '/'):
        if 'runinv' in simname:
            
            ctype = thiscluster.split("_")[1] #get all, FoV, etc
            filename = (plotconfig.outpath + '/' + simname + '/' +
                        thiscluster + '/c_of_m_3D.dat')
            
            macro = np.loadtxt(filename)
            nsnap =  macro[:,0]
            lagrange = macro[:,4]
            time = (nsnap/nsnap[-1])*duration
            #print ("   Doing ", simname, "...", sep="")
            plt.plot(time, lagrange)
            
    saveplot = (plotconfig.outpath + '/plots/Rh_'+ctype+'_comparek.pdf')
    plt.tight_layout()
    plt.savefig(saveplot, bbox_inches='tight')
    print ('         Rh_'+ctype+'_comparek.pdf')
    print ("")
    plt.close()
