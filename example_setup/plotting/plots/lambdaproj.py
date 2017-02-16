#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container #needed to remove error bars from legend
from matplotlib.ticker import AutoMinorLocator #for minor ticks
from sys import argv

fbin, fdim_val, qvir_val, outpath = argv[1:5]

path = outpath + '/outputs/'

duration = 10 #duration of simulation (Myr)


for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        for clustype in os.listdir(path + '/' + simname + '/'):
            if 'cluster_' in clustype:

                kval = simname.split("_")[1] #get k01, k02, etc
                clusterstring = clustype.split("_")[1] #'all', 'FoV#', 'rhalf#'
                
                filepath = (path + simname + '/' + clustype + '/lambda')
                
                #Get each type of lambda:
                #Use file names for xy projections (most likely projection to have calculated)
                for lamfile in os.listdir(filepath + '/'):
                    if '_xy' in lamfile:
                        
                        lamtype = lamfile.split("_")[1] #lambar, lamrms, etc
                        
                        #Find file names in xy, yz, xz, 3D:
                        filenamexy = (filepath + '/MST_' + lamtype + '_xy.dat')
                        filenameyz = (filepath + '/MST_' + lamtype + '_yz.dat')
                        filenamexz = (filepath + '/MST_' + lamtype + '_xz.dat')
                        filename3D = (filepath + '/MST_' + lamtype + '_3D.dat')
                        #Load the data from the text files:
                        lamdat_xy = np.loadtxt(filenamexy)
                        lamdat_yz = np.loadtxt(filenameyz)
                        lamdat_xz = np.loadtxt(filenamexz)
                        lamdat_3D = np.loadtxt(filename3D)
                        
                        #Find the time of each snapshot:
                        #(Only needs to be done once but easier to keep in loop as need filename)
                        nsnap =  lamdat_xy[:,0]
                        time = (nsnap/nsnap[-1])*duration
                        
                        #Save lambda data & errors for each projection in lists:
                        lam_xy = lamdat_xy[:,3]
                        lam_low_xy = lamdat_xy[:,4]
                        lam_up_xy = lamdat_xy[:,5]
                        yerr_low_xy = lam_xy - lam_low_xy
                        yerr_up_xy = lam_up_xy - lam_xy
                        
                        lam_yz = lamdat_yz[:,3]
                        lam_low_yz = lamdat_yz[:,4]
                        lam_up_yz = lamdat_yz[:,5]
                        yerr_low_yz = lam_yz - lam_low_yz
                        yerr_up_yz = lam_up_yz - lam_yz
                        
                        lam_xz = lamdat_xz[:,3]
                        lam_low_xz = lamdat_xz[:,4]
                        lam_up_xz = lamdat_xz[:,5]
                        yerr_low_xz = lam_xz - lam_low_xz
                        yerr_up_xz = lam_up_xz - lam_xz
                        
                        lam_3D = lamdat_3D[:,3]
                        lam_low_3D = lamdat_3D[:,4]
                        lam_up_3D = lamdat_3D[:,5]
                        yerr_low_3D = lam_3D - lam_low_3D
                        yerr_up_3D = lam_up_3D - lam_3D
                        
                        
                        #==================================================
                        #                   Make the plots!
                        #==================================================
                        
                        my_dpi=96
                        #set minor tick marks:
                        minorLocator = AutoMinorLocator(4) #this gives n-1 minor ticks
                        #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
                        
                        plt.figure()
                        plt.xlabel("Time (Myr)")
                        if "lambar" in lamtype:
                            plt.ylabel(r"$\overline{\Lambda}$")
                        elif "lamrms" in lamtype:
                            plt.ylabel(r"$\Lambda_{rms}$")
                        elif "lamsmr" in lamtype:
                            plt.ylabel(r"$\Lambda_{smr}$")
                        elif "lamhar" in lamtype:
                            plt.ylabel(r"$\Lambda_{har}$")
                        elif "lamtil" in lamtype:
                            plt.ylabel(r"$\widetilde{\Lambda}$")
                        elif "lamstar" in lamtype:
                            plt.ylabel(r"$\Lambda^\star}$")
                        elif "gam" in lamtype:
                            plt.ylabel(r"$\Gamma$")
                        elif "lamln" in lamtype:
                            plt.ylabel(r"$ln(\Lambda)$")
                        else:
                            plt.ylabel(lamtype)
                            #('r' flag is 'raw string literal')
                            
                        plt.title("Mass segregation")
                        #plt.plot(time, lam_xy)
                        plt.errorbar(time,lam_xy,yerr=[yerr_low_xy,yerr_up_xy],errorevery=5,label='xy')
                        plt.errorbar(time,lam_yz,yerr=[yerr_low_yz,yerr_up_yz],errorevery=5,label='yz')
                        plt.errorbar(time,lam_xz,yerr=[yerr_low_xz,yerr_up_xz],errorevery=5,label='xz')
                        plt.errorbar(time,lam_3D,yerr=[yerr_low_3D,yerr_up_3D],errorevery=5,label='3D')
                        plt.ylim(0,15)
                        plt.axes().xaxis.set_minor_locator(minorLocator)
                        
                        plt.annotate("fdim = " + str(fdim_val) + ", qvir = " + str(qvir_val), 
                                     xy=(0.98, 1.01), xycoords='axes fraction', 
                                     horizontalalignment='right', verticalalignment='bottom', 
                                     fontsize=10)
                        
                        plt.annotate("fbin = " + str(fbin), xy=(0., -0.09),
                                     xycoords='axes fraction', horizontalalignment='right',
                                     verticalalignment='bottom', fontsize=10)
                        
                        plt.annotate(str(kval), xy=(0.99, -0.08),
                                     xycoords='axes fraction', horizontalalignment='right',
                                     verticalalignment='bottom', fontsize=10)
                        
                        plt.legend(loc="upper right", fontsize=10)
                        
                        saveplot = (path + 'plots/' + kval + '_' + 
                                    clusterstring + '_' + lamtype + '.pdf')
                        plt.tight_layout()
                        plt.savefig(saveplot, bbox_inches='tight')
                        
                        print "Graph saved at " + saveplot
                        plt.close()
