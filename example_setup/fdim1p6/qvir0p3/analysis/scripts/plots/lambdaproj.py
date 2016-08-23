#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container #needed to remove error bars from legend
from sys import argv

fbin, fdim_val, qvir_val = argv[1:4]

path = '../outputs'

duration = 10 #duration of simulation (Myr)


for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        for clustype in os.listdir(path + '/' + simname + '/'):
            if 'cluster_' in clustype:

                kval = simname.split("_")[1] #get k01, k02, etc

                filenamexy = path + '/' + simname + '/' + clustype + '/lambda_xy'
                filenameyz = path + '/' + simname + '/' + clustype + '/lambda_yz'
                filenamexz = path + '/' + simname + '/' + clustype + '/lambda_xz'
                filename3D = path + '/' + simname + '/' + clustype + '/lambda_3D'

                clusterstring = clustype.split("_")[1]

                lambdxy = np.loadtxt(filenamexy)
                lambdyz = np.loadtxt(filenameyz)
                lambdxz = np.loadtxt(filenamexz)
                lambd3D = np.loadtxt(filename3D)

                nsnap =  lambd3D[:,0]
                time = (nsnap/nsnap[-1])*duration

                #Lambda bar:

                lam_barxy = lambdxy[:,1]
                lam_low_barxy = lambdxy[:,2]
                lam_up_barxy = lambdxy[:,3]
                yerr_low_barxy = lam_barxy - lam_low_barxy
                yerr_up_barxy = lam_up_barxy - lam_barxy

                lam_baryz = lambdyz[:,1]
                lam_low_baryz = lambdyz[:,2]
                lam_up_baryz = lambdyz[:,3]
                yerr_low_baryz = lam_baryz - lam_low_baryz
                yerr_up_baryz = lam_up_baryz - lam_baryz

                lam_barxz = lambdxz[:,1]
                lam_low_barxz = lambdxz[:,2]
                lam_up_barxz = lambdxz[:,3]
                yerr_low_barxz = lam_barxz - lam_low_barxz
                yerr_up_barxz = lam_up_barxz - lam_barxz

                lam_bar3D = lambd3D[:,1]
                lam_low_bar3D = lambd3D[:,2]
                lam_up_bar3D = lambd3D[:,3]
                yerr_low_bar3D = lam_bar3D - lam_low_bar3D
                yerr_up_bar3D = lam_up_bar3D - lam_bar3D

                #Lambda tilde:

                lam_tilxy = lambdxy[:,4]
                lam_low_tilxy = lambdxy[:,5]
                lam_up_tilxy = lambdxy[:,6]
                yerr_low_tilxy = lam_tilxy - lam_low_tilxy
                yerr_up_tilxy = lam_up_tilxy - lam_tilxy

                lam_tilyz = lambdyz[:,4]
                lam_low_tilyz = lambdyz[:,5]
                lam_up_tilyz = lambdyz[:,6]
                yerr_low_tilyz = lam_tilyz - lam_low_tilyz
                yerr_up_tilyz = lam_up_tilyz - lam_tilyz

                lam_tilxz = lambdxz[:,4]
                lam_low_tilxz = lambdxz[:,5]
                lam_up_tilxz = lambdxz[:,6]
                yerr_low_tilxz = lam_tilxz - lam_low_tilxz
                yerr_up_tilxz = lam_up_tilxz - lam_tilxz

                lam_til3D = lambd3D[:,4]
                lam_low_til3D = lambd3D[:,5]
                lam_up_til3D = lambd3D[:,6]
                yerr_low_til3D = lam_til3D - lam_low_til3D
                yerr_up_til3D = lam_up_til3D - lam_til3D

                #Lambda star:

                lam_starxy = lambdxy[:,7]
                lam_low_starxy = lambdxy[:,8]
                lam_up_starxy = lambdxy[:,9]
                yerr_low_starxy = lam_starxy - lam_low_starxy
                yerr_up_starxy = lam_up_starxy - lam_starxy

                lam_staryz = lambdyz[:,7]
                lam_low_staryz = lambdyz[:,8]
                lam_up_staryz = lambdyz[:,9]
                yerr_low_staryz = lam_staryz - lam_low_staryz
                yerr_up_staryz = lam_up_staryz - lam_staryz

                lam_starxz = lambdxz[:,7]
                lam_low_starxz = lambdxz[:,8]
                lam_up_starxz = lambdxz[:,9]
                yerr_low_starxz = lam_starxz - lam_low_starxz
                yerr_up_starxz = lam_up_starxz - lam_starxz

                lam_star3D = lambd3D[:,7]
                lam_low_star3D = lambd3D[:,8]
                lam_up_star3D = lambd3D[:,9]
                yerr_low_star3D = lam_star3D - lam_low_star3D
                yerr_up_star3D = lam_up_star3D - lam_star3D

                #Gamma:

                gamxy = lambdxy[:,10]
                gam_lowxy = lambdxy[:,11]
                gam_upxy = lambdxy[:,12]
                yerr_low_gamxy = gamxy - gam_lowxy
                yerr_up_gamxy = gam_upxy - gamxy

                gamyz = lambdyz[:,10]
                gam_lowyz = lambdyz[:,11]
                gam_upyz = lambdyz[:,12]
                yerr_low_gamyz = gamyz - gam_lowyz
                yerr_up_gamyz = gam_upyz - gamyz

                gamxz = lambdxz[:,10]
                gam_lowxz = lambdxz[:,11]
                gam_upxz = lambdxz[:,12]
                yerr_low_gamxz = gamxz - gam_lowxz
                yerr_up_gamxz = gam_upxz - gamxz

                gam3D = lambd3D[:,10]
                gam_low3D = lambd3D[:,11]
                gam_up3D = lambd3D[:,12]
                yerr_low_gam3D = gam3D - gam_low3D
                yerr_up_gam3D = gam_up3D - gam3D


                my_dpi=96
                #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
                
                #==================================================
                #                   Make the plots!
                #==================================================
                #
                #------------#
                # Lambda bar #
                #------------#

                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(r"$\overline{\Lambda}$")
                plt.title("Mass segregation")
                #plt.plot(time, lam_barxy)
                plt.errorbar(time,lam_barxy,yerr=[yerr_low_barxy,yerr_up_barxy],errorevery=5,label='xy')
                plt.errorbar(time,lam_baryz,yerr=[yerr_low_baryz,yerr_up_baryz],errorevery=5,label='yz')
                plt.errorbar(time,lam_barxz,yerr=[yerr_low_barxz,yerr_up_barxz],errorevery=5,label='xz')
                plt.errorbar(time,lam_bar3D,yerr=[yerr_low_bar3D,yerr_up_bar3D],errorevery=5,label='3D')
                plt.ylim(0,15)

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

                saveplot = (path + '/plots/' + kval + '_' + 
                            clusterstring + '_lambar.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                print "Graph saved at " + saveplot
                plt.close()


                #--------------#
                # Lambda tilde #
                #--------------#

                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(r"$\widetilde{\Lambda}$")
                plt.title("Mass segregation")
                #plt.plot(time, lam_til)
                plt.errorbar(time,lam_tilxy,yerr=[yerr_low_tilxy,yerr_up_tilxy],errorevery=5,label='xy')
                plt.errorbar(time,lam_tilyz,yerr=[yerr_low_tilyz,yerr_up_tilyz],errorevery=5,label='yz')
                plt.errorbar(time,lam_tilxz,yerr=[yerr_low_tilxz,yerr_up_tilxz],errorevery=5,label='xz')
                plt.errorbar(time,lam_til3D,yerr=[yerr_low_til3D,yerr_up_til3D],errorevery=5,label='3D')
                plt.ylim(0,25)

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

                saveplot = (path + '/plots/' + kval + '_' + 
                            clusterstring + '_lamtil.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                print "Graph saved at " + saveplot
                plt.close()


                #-------------#
                # Lambda star #
                #-------------#

                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(r"$\Lambda^\star$")
                plt.title("Mass segregation")
                #plt.plot(time, lam_star)
                plt.errorbar(time,lam_starxy,yerr=[yerr_low_starxy,yerr_up_starxy],errorevery=5,label='xy')
                plt.errorbar(time,lam_staryz,yerr=[yerr_low_staryz,yerr_up_staryz],errorevery=5,label='yz')
                plt.errorbar(time,lam_starxz,yerr=[yerr_low_starxz,yerr_up_starxz],errorevery=5,label='xz')
                plt.errorbar(time,lam_star3D,yerr=[yerr_low_star3D,yerr_up_star3D],errorevery=5,label='3D')
                plt.ylim(0,15)

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

                saveplot = (path + '/plots/' + kval + '_' + 
                            clusterstring + '_lamstar.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                print "Graph saved at " + saveplot
                plt.close()
                
                
                #-----------#
                #   Gamma   #
                #-----------#

                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(r"$\Gamma$")
                plt.title("Mass segregation")
                #plt.plot(time, gam)
                plt.errorbar(time,gamxy,yerr=[yerr_low_gamxy,yerr_up_gamxy],errorevery=5,label='xy')
                plt.errorbar(time,gamyz,yerr=[yerr_low_gamyz,yerr_up_gamyz],errorevery=5,label='yz')
                plt.errorbar(time,gamxz,yerr=[yerr_low_gamxz,yerr_up_gamxz],errorevery=5,label='xz')
                plt.errorbar(time,gam3D,yerr=[yerr_low_gam3D,yerr_up_gam3D],errorevery=5,label='3D')
                plt.ylim(0,25)

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

                saveplot = (path + '/plots/' + kval + '_' + 
                            clusterstring + '_gamma.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                print "Graph saved at " + saveplot
                plt.close()
