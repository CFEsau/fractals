#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

#print len(argv)
fbin, fdim_val, qvir_val = argv[1:4]

#path = str(raw_input('input the path to the simulation directory: '))
path = '../outputs'
#print path

#fdim_val = raw_input("fdim: ")
#qvir_val = raw_input("qvir: ")
#duration = float(input('input the duration of the simulation (Myr): '))
duration = 10


for simname in os.listdir(path + '/'):
    if 'runinv' in simname:

        kval = simname.split("_")[1] #get k01, k02, etc
        filename = path + '/' + simname + '/lambda'

        lambd = np.loadtxt(filename)
        nsnap =  lambd[:,0]
        time = (nsnap/nsnap[-1])*duration
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
