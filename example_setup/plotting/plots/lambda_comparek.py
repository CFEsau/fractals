#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

fbin, fdim_val, qvir_val, outpath = argv[1:5] #use when defining parameters inline

path = outpath + '/outputs/'

duration = 10. #Duration of simulation (Myr)

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

#find a k value that is present in outputs directory
#(use idir instead of os.listdir(path)[0] as 'plots' directory comes first)
idir=0
for alldir in os.listdir(path):
    if 'runinv_k' in alldir:
        firstk=os.listdir(path)[idir]
        #print ("first k directory: ",firstk)
        break
    idir=idir+1

#make list for names of different lambda:
lambda_types=[]
lambda_types.append('lambar')
lambda_types.append('lamrms')
lambda_types.append('lamsmr')
lambda_types.append('lamhar')
lambda_types.append('lamtil')
lambda_types.append('lamstar')
lambda_types.append('gam')
lambda_types.append('lamln')

#dictionary for Y/N logical of what to plot:
lambda_to_plot=dict()
lambda_to_plot['lambar']= 'Y'
lambda_to_plot['lamrms']= 'Y'
lambda_to_plot['lamsmr']= 'N'
lambda_to_plot['lamhar']= 'N'
lambda_to_plot['lamtil']= 'Y'
lambda_to_plot['lamstar']='N'
lambda_to_plot['gam']=    'Y'
lambda_to_plot['lamln']=  'Y'

#dictionary for tex strings:
lambda_latex=dict()
lambda_latex['lambar']='$\overline{\Lambda}$'
lambda_latex['lamrms']='$\Lambda_{rms}$'
lambda_latex['lamsmr']='$\Lambda_{smr}$'
lambda_latex['lamhar']='$\Lambda_{har}$'
lambda_latex['lamtil']='$\widetilde{\Lambda}$'
lambda_latex['lamstar']='$\Lambda^\star$'
lambda_latex['gam']='$\Gamma$'
lambda_latex['lamln']='$\mathrm{ln}(\Lambda)$'


for clustertype in os.listdir(path + '/' + firstk + '/'):
    if 'cluster' in clustertype:
        for lamtype in lambda_types:
            if lambda_to_plot[lamtype] == 'Y':
                print ("   Doing",lamtype,"...", sep=" ")
                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(r"$\overline{\Lambda}$")
                plt.title("Mass segregation")
                plt.ylim(0,15)
                #plt.xticks(np.arange(min(x), max(x), 1.0))
                #plt.locator_params(axis='x',nbins=10)
                plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                             xycoords='axes fraction', fontsize = 10)
                plt.annotate('fdim = ' + fdim_val + ', qvir = ' + qvir_val, 
                             xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
                
                for simname in os.listdir(path + '/'):
                    if 'runinv' in simname:
                        
                        filename = (path + '/' + simname + '/' +
                                    clustertype + '/lambda/MST_'+lamtype+'_3D.dat')
                        cstring = clustertype.split("_")[1] #get all, FoV, etc
                        
                        lambd = np.loadtxt(filename)
                        nsnap =  lambd[:,0]
                        time = (nsnap/nsnap[-1])*duration
                        lamdata = lambd[:,3]
                        
                        print ("       Plotting ", simname, "...", sep="")
                        plt.plot(time, lamdata)
                        
                saveplot = path + 'plots/'+lamtype+'_'+cstring+'_allk.pdf'
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                print ("Saved in ", saveplot)
                plt.close()
