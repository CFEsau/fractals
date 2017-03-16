#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator #for minor ticks
from sys import argv

fbin, fdim_val, qvir_val, outpath = argv[1:5] #use when defining parameters inline

path = outpath + '/outputs/'

duration = 10. #duration of simulation (Myr)

my_dpi=96
#plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 sheet is 8.3 x 11.7 (portrait)
saveplot = ''

#set minor tick marks:
minorLocator = AutoMinorLocator(4) #this gives n-1 minor ticks

#===================================================
#  Comparing different measures of mass segregation
#===================================================


for simname in os.listdir(path + '/'):
    if 'runinv' in simname: #runinv_k01, etc
        for clustype in os.listdir(path + '/' + simname + '/'): #all, FoV, etc
            if 'cluster_' in clustype:
                
                kval = simname.split("_")[1] #get k01, k02, etc
                
                #set up plot:
                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(r"$\Lambda$")
                plt.title("Mass segregation")
                plt.ylim(0,25)
                #plt.xticks(np.arange(0, 11, 1.0))
                #plt.locator_params(axis='x',nbins=10)
                plt.axes().xaxis.set_minor_locator(minorLocator)

                #These next bits are to plot different lambda in certain order:
                
                #make list for names of different lambda:
                lambda_types=[]
                lambda_types.append('lambar') #arithmetic mean
                lambda_types.append('lamrms') #rms
                lambda_types.append('lamsmr') #smr
                lambda_types.append('lamhar') #harmonic mean
                lambda_types.append('lamtil') #median
                lambda_types.append('lamNmed') #mean of N median (~2 or 3)
                lambda_types.append('lamstar') #weird one... ignore
                lambda_types.append('gam') #geometric mean
                lambda_types.append('lamln') #log

                #dictionary for Y/N logical of what to plot:
                lambda_to_plot=dict()
                lambda_to_plot['lambar']= 'Y'
                lambda_to_plot['lamrms']= 'Y'
                lambda_to_plot['lamsmr']= 'N'
                lambda_to_plot['lamhar']= 'N'
                lambda_to_plot['lamtil']= 'Y'
                lambda_to_plot['lamNmed']='N'
                lambda_to_plot['lamstar']='N'
                lambda_to_plot['gam']=    'Y'
                lambda_to_plot['lamln']=  'N'

                #dictionary for lambda data:
                lambda_data=dict()

                #dictionary for tex strings:
                lambda_tex=dict()
                lambda_tex['lambar']='$\overline{\Lambda}$'
                lambda_tex['lamrms']='$\Lambda_{rms}$'
                lambda_tex['lamsmr']='$\Lambda_{smr}$'
                lambda_tex['lamhar']='$\Lambda_{har}$'
                lambda_tex['lamtil']='$\widetilde{\Lambda}$'
                lambda_tex['lamNmed']='$\widetilde{\Lambda}_N$'
                lambda_tex['lamstar']='$\Lambda^\star$'
                lambda_tex['gam']='$\Gamma$'
                lambda_tex['lamln']='$\mathrm{ln}(\Lambda)$'
                
                #load lambda data:
                lampath = path + '/' + simname + '/' + clustype + '/lambda/'
                
                for lamfile in os.listdir(lampath):
                    if '_3D' in lamfile:
                        lamtype = lamfile.split("_")[1] #get lambda type (bar, rms, etc)
                        lam = np.loadtxt(lampath + lamfile) #load data
                        nsnap =  lam[:,0]
                        time = (nsnap/nsnap[-1])*duration
                        lambda_data[lamtype] = lam[:,3] # save data for each lamtype
                
                for lamtype in lambda_types:
                    if lambda_to_plot[lamtype] == 'Y':
                        plt.plot(time,lambda_data[lamtype],label=lambda_tex[lamtype])
                        plt.legend(fontsize=10)
                        plt.annotate("fbin = " + fbin, xy=(0.95,-0.09), 
                                     xycoords='axes fraction', fontsize = 10)
                        plt.annotate('fdim = ' + fdim_val + ', qvir = ' + qvir_val, 
                                     xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
                        
                        cstring = clustype.split("_")[1]
                    elif lambda_to_plot[lamtype] == 'N':
                        print (" ",lamtype," not plotted")
                    else:
                        print ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                        print ("WARNING: Invalid plot option '" + lambda_to_plot[ltype] + "' entered !" + "for" + lamtype)
                
                #plt.show()
                saveplot = (path + 'plots/' + kval + '_' + cstring + '_alllam.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                print ("   Graph saved at ", saveplot)
                plt.close()
