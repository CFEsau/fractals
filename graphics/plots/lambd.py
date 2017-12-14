#!/usr/bin/env python

#from __future__ import print_function
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import container #needed to remove error bars from legend
from matplotlib.ticker import AutoMinorLocator #for minor ticks
import plotconfig
    
def printlambda():
    print "\n------------"
    print "-- Lambda --"
    print "------------"

def lambdasetup():
    #This plots lambda values in user-specified
    #order, not default alphabetically.
    #make list for names of different lambda:
    plotconfig.lambdatypes=[]
    plotconfig.lambdatypes.append('lambar') #arithmetic mean
    plotconfig.lambdatypes.append('lamrms') #rms
    plotconfig.lambdatypes.append('lamsmr') #smr
    plotconfig.lambdatypes.append('lamhar') #harmonic mean
    plotconfig.lambdatypes.append('lamtil') #median
    plotconfig.lambdatypes.append('lam3med') #mean of N median (~2 or 3)
    plotconfig.lambdatypes.append('lamstar') #weird one... ignore
    plotconfig.lambdatypes.append('gam') #geometric mean
    plotconfig.lambdatypes.append('lamln') #log

    #dictionary of what to plot (Y/N):
    lambda_to_plot={}
    lambda_to_plot['lambar']= 'Y'
    lambda_to_plot['lamrms']= 'Y'
    lambda_to_plot['lamsmr']= 'N'
    lambda_to_plot['lamhar']= 'N'
    lambda_to_plot['lamtil']= 'N'
    lambda_to_plot['lam3med']='Y'
    lambda_to_plot['lamstar']='N'
    lambda_to_plot['gam']=    'Y'
    lambda_to_plot['lamln']=  'N'

    #update lambdatypes list to include only types with 'Y'
    plotconfig.lambdatypes = ([thislam for thislam in plotconfig.lambdatypes
                               if ('Y' in lambda_to_plot[thislam])])
    for thislam in plotconfig.lambdatypes:
        print (thislam,lambda_to_plot[thislam])
        
    #dictionary for tex strings:
    plotconfig.lambda_tex={}
    plotconfig.lambda_tex['lambar']='$\overline{\Lambda}$'
    plotconfig.lambda_tex['lamrms']='$\Lambda_{rms}$'
    plotconfig.lambda_tex['lamsmr']='$\Lambda_{smr}$'
    plotconfig.lambda_tex['lamhar']='$\Lambda_{har}$'
    plotconfig.lambda_tex['lamtil']='$\widetilde{\Lambda}$'
    plotconfig.lambda_tex['lam3med']='$\widetilde{\Lambda}_3$'
    plotconfig.lambda_tex['lamstar']='$\Lambda^\star$'
    plotconfig.lambda_tex['gam']='$\Gamma$'
    plotconfig.lambda_tex['lamln']='$\mathrm{ln}(\Lambda)$'
    

def lambdaprojections(thiscluster):
    #duration = 10. #Duration of simulation (Myr)
    projections=['3D','xy','xz','yz']
    
    for simname in os.listdir(plotconfig.outpath + '/'):
        #loop through each simulation (k number):
        if 'runinv' in simname:
            kval = simname.split("_")[1] #get k01, k02, etc
            ctype = thiscluster.split("_")[1] #'all', FoV, etc
            filepath = (plotconfig.outpath+'/'+simname+'/' + 
                        thiscluster + '/lambda')
            
            for thislambda in plotconfig.lambdatypes:
                #print("Doing",thislambda)
                #check that this file exists:
                lamfilelist=[]
                found=''
                for filename in os.listdir(filepath):
                    if thislambda in filename:
                        #print thislambda
                        found = 'found'
                        lamfilelist.append(filename)
                        #break #had this before moving 'lamfilelist' to here        
                if not found:
                    print "WARNING: %s doesn't exist. Skipping." % thislambda
                
                #read in data & make plots for this lambda
                #set up plot:
                saveplot = ''
                my_dpi = 96
                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(plotconfig.lambda_tex[thislambda])
                plt.title("Mass segregation")
                #n-1 minor ticks:
                plt.axes().xaxis.set_minor_locator(AutoMinorLocator(4))
                plt.ylim(0,15)
                plt.annotate("fdim = "+plotconfig.fdim+", qvir = "+plotconfig.qvir,
                             xy=(1.0, -0.09), xycoords='axes fraction',
                             horizontalalignment='right',
                             verticalalignment='bottom', fontsize=10)
                plt.annotate("fbin = " + plotconfig.fbin, xy=(0.05, -0.09),
                             xycoords='axes fraction', horizontalalignment='right',
                             verticalalignment='bottom', fontsize=10)
                plt.annotate(kval, xy=(0.99, 1.01),
                             xycoords='axes fraction', horizontalalignment='right',
                             verticalalignment='bottom', fontsize=10)
                
                #plt.show()
                #do this if having the 'if errorbars' inside
                #the 'do' takes ages
                #if plotconfig.errorbars=='y':
                #    makeplot='plt.errorbar(time,lambda_data, \
                #    yerr=[lambda_minerr,lambda_maxerr], \
                #    errorevery=5,label=thisproj)'
                #else:
                #    makeplot='plt.plot(time,lambda_data,label=thisproj)'
                
                #plot data for each projection:
                for thisfile in lamfilelist:
                    for thisproj in projections:
                        #print thisproj
                        if thisproj in thisfile:
                            #for thisfile in lamfilelist:
                            nsnap = np.loadtxt(filepath+'/'+thisfile)[:,0]
                            time = (nsnap/nsnap[-1])*plotconfig.duration
                            
                            lambda_data = np.loadtxt(filepath+'/'+thisfile)[:,3]
                            #if thislambda == 'lambar' and kval=='k01':
                            #    print thisfile
                            #    print lambda_data[0:3]
                            lambda_min = np.loadtxt(filepath+'/'+thisfile)[:,4]
                            lambda_max = np.loadtxt(filepath+'/'+thisfile)[:,5]
                            lambda_maxerr = lambda_data - lambda_min
                            lambda_minerr = lambda_max - lambda_data
                            #make plot
                            #exec(makeplot) #if 'if errorbars' loop takes ages
                            if plotconfig.errorbars=='y':
                                plt.errorbar(time,lambda_data,
                                             yerr=[lambda_minerr,lambda_maxerr],
                                             errorevery=5,label=thisproj)
                            else:
                                plt.plot(time,lambda_data,label=thisproj)
                            break #skip next projections & move on to next file
                        
                plt.legend(loc="upper right", fontsize=10)
                #plt.show()
                saveplot = (plotconfig.outpath+'/plots/'+
                                kval+'_'+ctype+'_'+thislambda+'.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                #print ("    ",kval,"_",ctype,"_",thislambda,".pdf",sep="")
                print "\t %s_%s_%s.pdf" % (kval, ctype, thislambda)
                plt.close()

def projectioncompare(thiscluster):
    #duration = 10. #Duration of simulation (Myr)
    projections2D=['xy','xz','yz']
    
    for simname in os.listdir(plotconfig.outpath + '/'):
        #loop through each simulation (k number):
        if 'runinv' in simname:
            kval = simname.split("_")[1] #get k01, k02, etc
            ctype = thiscluster.split("_")[1] #'all', FoV, etc
            filepath = (plotconfig.outpath+'/'+simname+'/' + 
                        thiscluster + '/lambda')
            
            for thislambda in plotconfig.lambdatypes:
                #print "Doing %s" % thislambda
                #check that this file exists:
                lamfilelist=[]
                found=''
                for filename in os.listdir(filepath):
                    if thislambda in filename:
                        #print thislambda
                        found = 'found'
                        lamfilelist.append(filename)
                        #break #had this before moving 'lamfilelist' to here        
                if not found:
                    print "WARNING: %s doesn't exist. Skipping." % thislambda
                
                #read in data & make plots for this lambda
                #set up plot:
                saveplot = ''
                my_dpi = 96
                plt.figure()
                plt.xlabel("Time (Myr)")
                plt.ylabel(plotconfig.lambda_tex[thislambda])
                plt.title("Mass segregation for 2D projections relative to 3D")
                #n-1 minor ticks:
                plt.axes().xaxis.set_minor_locator(AutoMinorLocator(4))
                plt.annotate("fdim = "+plotconfig.fdim+", qvir = "+plotconfig.qvir,
                             xy=(1.0, -0.09), xycoords='axes fraction',
                             horizontalalignment='right',
                             verticalalignment='bottom', fontsize=10)
                plt.annotate("fbin = " + plotconfig.fbin, xy=(0.05, -0.09),
                             xycoords='axes fraction', horizontalalignment='right',
                             verticalalignment='bottom', fontsize=10)
                plt.annotate(kval, xy=(0.99, 1.01),
                             xycoords='axes fraction', horizontalalignment='right',
                             verticalalignment='bottom', fontsize=10)
                
                ylow=0
                #plot data for each 2D projection:
                #'thisfile' is MST_[thislambda]_[proj].dat
                for thisfile in lamfilelist:
                    #read in 3D data:
                    if '3D' in thisfile:
                        nsnap = np.loadtxt(filepath+'/'+thisfile)[:,0]
                        time = (nsnap/nsnap[-1])*plotconfig.duration
                            
                        data_3D = np.loadtxt(filepath+'/'+thisfile)[:,3]
                        lambda_min = np.loadtxt(filepath+'/'+thisfile)[:,4]
                        lambda_max = np.loadtxt(filepath+'/'+thisfile)[:,5]
                        lambda_maxerr = data_3D - lambda_min
                        lambda_minerr = lambda_max - data_3D
                        #make plot
                        zero_line_data = np.array([0 for i in xrange(len(time))])
                        plt.plot(time,data_3D,label="3D MST",linestyle='dotted')
                        plt.plot(time,zero_line_data)
                        
                    for thisproj in projections2D:
                        #print thisproj
                        if thisproj in thisfile:
                            #for thisfile in lamfilelist:
                            nsnap = np.loadtxt(filepath+'/'+thisfile)[:,0]
                            time = (nsnap/nsnap[-1])*plotconfig.duration
                            
                            lambda_data = np.loadtxt(filepath+'/'+thisfile)[:,3]
                            #if thislambda == 'lambar' and kval=='k01':
                            #    print thisfile
                            #    print lambda_data[0:3]
                            lambda_min = np.loadtxt(filepath+'/'+thisfile)[:,4]
                            lambda_max = np.loadtxt(filepath+'/'+thisfile)[:,5]
                            lambda_maxerr = lambda_data - lambda_min
                            lambda_minerr = lambda_max - lambda_data
                            #make plot
                            plt.plot(time,lambda_data-data_3D,label=thisproj+'-3D')
                            ylow=min(min(lambda_data-data_3D),ylow)
                            #break #skip next projections & move on to next file

                #find minimum y value for axis limit:
                plt.ylim(ylow,15)
                plt.legend(loc="upper right", fontsize=10)
                #plt.show()
                saveplot = (plotconfig.outpath+'/plots/'+
                                kval+'_'+ctype+'_'+thislambda+'proj.pdf')
                plt.tight_layout()
                plt.savefig(saveplot, bbox_inches='tight')
                #print ("    ",kval,"_",ctype,"_",thislambda,"proj.pdf",sep="")
                print ("\t %s_%s_%sproj.pdf" % (kval, ctype, thislambda))
                plt.close()

#def 3Ddivergence(thiscluster):
    #find mass segregation divergence and dt divergence
    #use lambda > 2
    #if |2D-3D|>2, bounded by t1 and t2. dt=t2-t1
    # XXXXX no as don't have value for lambda_div. can't just use max.
    #use area under curve?
    # & how many points should we take for the scatter plot?
    

def lambdacompare(thiscluster):
    for simname in os.listdir(plotconfig.outpath + '/'):
        #loop through each simulation (k number):
        if 'runinv' in simname:
            kval = simname.split("_")[1] #get k01, k02, etc
            ctype = thiscluster.split("_")[1] #'all', FoV, etc
            filepath = (plotconfig.outpath+'/'+simname+'/' + 
                        thiscluster + '/lambda')
            
            #set up plot:
            saveplot = ''
            my_dpi = 96
            minorLocator = AutoMinorLocator(4) #n-1 minor tick marks
            plt.figure()
            plt.xlabel("Time (Myr)")
            plt.ylabel(r"$\Lambda$")
            plt.title("Different measures of mass segregation")
            #n-1 minor ticks:
            plt.axes().xaxis.set_minor_locator(minorLocator)
            plt.ylim(0,15)
            plt.annotate("fdim = "+plotconfig.fdim+", qvir = "+plotconfig.qvir,
                         xy=(1.0, -0.09), xycoords='axes fraction',
                         horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate("fbin = " + plotconfig.fbin, xy=(0.05, -0.09),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            plt.annotate(kval, xy=(0.99, 1.01),
                         xycoords='axes fraction', horizontalalignment='right',
                         verticalalignment='bottom', fontsize=10)
            #plt.show() #check setup

            
            #thislambda is given in lambdatypes
            #lamfilelist=[]
            for thislambda in plotconfig.lambdatypes:
                #Make sure file exists file exists:
                found=''
                for thisfile in glob.glob(filepath+'/*3D*'):
                    #make list of files to plot:
                    #(need list to make sure they're in order)
                    if thislambda in thisfile:
                        found = 'found'
                        #lamfilelist.append(thisfile)
                        
                        #read in data and plot:
                        nsnap = np.loadtxt(thisfile)[:,0]
                        time = (nsnap/nsnap[-1])*plotconfig.duration
                        lambda_data = np.loadtxt(thisfile)[:,3]
                        lambda_min = np.loadtxt(thisfile)[:,4]
                        lambda_max = np.loadtxt(thisfile)[:,5]
                        lambda_maxerr = lambda_data - lambda_min
                        lambda_minerr = lambda_max - lambda_data
                        
                        if plotconfig.errorbars=='y':
                            plt.errorbar(time,lambda_data,
                                 yerr=[lambda_minerr,lambda_maxerr],
                                 errorevery=5,
                                 label=plotconfig.lambda_tex[thislambda])
                        else:
                            plt.plot(time,lambda_data,
                                     label=plotconfig.lambda_tex[thislambda])
                        break #move on to next file
                    
                if not found:
                    print ("WARNING: %s data for this projection don't exist.",
                           " Skipping." % thislambda)
            
            plt.legend(loc="upper right", fontsize=10)
            #plt.show()
            saveplot = (plotconfig.outpath+'/plots/'+
                            kval+'_'+ctype+'_alllam.pdf')
            plt.tight_layout()
            plt.savefig(saveplot, bbox_inches='tight')
            #print ("    ",kval,"_",ctype,"_alllam.pdf",sep="")
            print "\t %s_%s_alllam.pdf" % (kval, ctype)
            plt.close()
            

def lambda_k(thiscluster):
    #Currently uses the 3D projection
    print "   Comparing across k..."
    
            #filepath = (plotconfig.outpath+'/'+simname+'/' + 
            #            thiscluster + '/lambda')
    ctype = thiscluster.split("_")[1] #get all, FoV, etc
    
    for thislambda in plotconfig.lambdatypes:
        #set up plot:
        my_dpi=96
        duration = 10. #Duration of simulation (Myr)
        #plt.figure(figsize=(11.7,8.3), dpi=my_dpi) # A4 landscape
        saveplot = ''
        plt.figure()
        plt.xlabel("Time (Myr)")
        plt.ylabel(plotconfig.lambda_tex[thislambda])
        plt.title("Mass segregation")
        plt.ylim(ymax = 15, ymin = 0)
        plt.annotate("fdim = "+plotconfig.fdim+
                     ", qvir = "+plotconfig.qvir,
                     xy=(1.0, -0.09),xycoords='axes fraction',
                     horizontalalignment='right',
                     verticalalignment='bottom', fontsize=10)
        plt.annotate("fbin = " + plotconfig.fbin, xy=(0.05, -0.09),
                     xycoords='axes fraction',
                     horizontalalignment='right',
                     verticalalignment='bottom', fontsize=10)
        
        for simname in os.listdir(plotconfig.outpath + '/'):
            if 'runinv' in simname:
                kval = simname.split("_")[1] #get k01, k02, etc
                filename = (plotconfig.outpath + '/' + simname + '/' +
                        thiscluster + '/lambda/MST_'+thislambda+'_3D.dat')
                #check file exists:
                if os.path.isfile(filename):
                    nsnap = np.loadtxt(filename)[:,0]
                    time = (nsnap/nsnap[-1])*plotconfig.duration
                    lambda_data = np.loadtxt(filename)[:,3]
                    lambda_min = np.loadtxt(filename)[:,4]
                    lambda_max = np.loadtxt(filename)[:,5]
                    lambda_maxerr = lambda_data - lambda_min
                    lambda_minerr = lambda_max - lambda_data
                    
                    print "\tPlotting %s..." % simname
                    plt.plot(time,lambda_data,label=kval)
                    #plt.show()
                    
                else:
                    print "WARNING: file for %s doesn't exist. Skipping." \
                        % (thislambda)
                    
        plt.legend(loc="upper right", fontsize=10)
        saveplot = (plotconfig.outpath + '/plots/'+thislambda+'_'+ctype+
                    '_comparek.pdf')
        #plt.show()
        plt.tight_layout()
        plt.savefig(saveplot, bbox_inches='tight')
        #print ("         ",thislambda,"_",ctype,"_comparek.pdf",sep="")
        print "\t %s_%s_comparek.pdf\n" % (thislambda, ctype)
        plt.close()
