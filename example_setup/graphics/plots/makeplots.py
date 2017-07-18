#!/usr/bin/env python

from __future__ import print_function
import os,sys,time
import plotconfig,filehandling
import energy,virial,rhalfm,lambd
#import timeit

skip='y'                   #skip energy, Q, rhalf plots? y/n
plotconfig.errorbars='n'   #plot error bars on lambda?
plotconfig.duration = 10   #Duration of simulation (Myr - for axis limits)

#Haven't made 'projection' a variable. Change manually. 3D at the moment.

plotconfig.fbin = ''
while not plotconfig.fbin:
    plotconfig.fbin = raw_input('Enter binary fraction: ')
    if plotconfig.fbin=='0%':
        fbindir='fbin0p0'
    elif plotconfig.fbin=='50%':
        fbindir='fbin0p5'
    elif plotconfig.fbin=='100%':
        fbindir='fbin1p0'
    else:
        print ("     Warning: Binary fraction",plotconfig.fbin,
               "not recognised. \n     Options are 0%, 50%, 100%.")
        plotconfig.fbin = ''

plotconfig.fdim = ''
while not plotconfig.fdim:
    plotconfig.fdim = raw_input('Enter fractal dimension: ')
    if plotconfig.fdim=='1.6':
        fdimdir='fdim1p6'
    elif plotconfig.fdim=='2.0':
        fdimdir='fdim2p0'
    elif plotconfig.fdim=='2.6':
        fdimdir='fdim2p6'
    elif plotconfig.fdim=='3.0':
        fdimdir='fdim3p0'
    else:
        print ("     Warning: Fractal dimension",plotconfig.fdim,
               "not recognised. \n     Options are 1.6, 2.0, 2.6, 3.0.")
        plotconfig.fdim = ''

plotconfig.qvir = ''
while not plotconfig.qvir:
    plotconfig.qvir = raw_input('Enter virial ratio: ')
    if plotconfig.qvir=='0.3':
        qvirdir='qvir0p3'
    elif plotconfig.qvir=='0.5':
        qvirdir='qvir0p5'
    else:
        print ("     Warning: Virial ratio",plotconfig.qvir,
               "not recognised. \n     Options are 0.3, 0.5.")
        plotconfig.qvir = ''

modelpath='../../'+fbindir+'/'+fdimdir+'/'+qvirdir
print ("")
print ("Model path is",modelpath)

print ("")
print ("**************************")
print ("*      Making plots      *")
print ("**************************")

plotconfig.outpath = modelpath + '/outputs'

#get different types of cluster (all, FoV, etc):
plotconfig.clustertypes = []
filehandling.getclusters()
#print ("Cluster list:",plotconfig.clustertypes)

#Make 'plots' directory if it doesn't exist:
if not os.path.exists(plotconfig.outpath+'/plots'):
    print ("Making new directory",plotconfig.outpath,"/plots...")
    print ("")
    os.makedirs(plotconfig.outpath+'/plots')

if skip == 'n':
    #----------------
    #     Energy
    #----------------
    energy.printenergy() #print banner
    energy.plotenergy() #plot energies for each model
    filehandling.mergefiles() #merge plots into one document
    energy.energy_k() #compare simulations (k## together)
    
    #----------------
    #      Qvir
    #----------------
    virial.printvirial() #print banner
    virial.plotvirial() #plot Q for each model
    filehandling.mergefiles() #merge plots into one document
    virial.virial_k() #compare simulations (k## together)
    
    #----------------
    #     Rhalf
    #----------------
    rhalfm.printrhalfm() #print banner
    #do for different cluster types:
    #print ("Cluster list:",plotconfig.clustertypes)
    for thiscluster in plotconfig.clustertypes:
        rhalfm.plotrhalfm(thiscluster) #plot rhalfm for each model
        filehandling.mergefiles() #merge plots into one document
        rhalfm.rhalfm_k(thiscluster) #compare simulations (k## together)
        print ("")
    
else:
    print("Skipping energy...")
    print("Skipping virial ratio...")
    print("Skipping half-mass radius...")

    
lambd.printlambda() #print banner
lambd.lambdasetup() #set which lambda types to plot & y-label text

#start_time = timeit.default_timer()
#overplot projections for all types of lambda, for each cluster

for thiscluster in plotconfig.clustertypes:
    #--------------------
    # Lambda projections
    #--------------------
    lambd.lambdaprojections(thiscluster) #compare projections for each lambda
    filehandling.mergefiles() #merge plots into one document
    print ("")
    
    #plot 2D projections relative to 3D for each lambda:
    lambd.projectioncompare(thiscluster)
    filehandling.mergefiles() #merge plots into one document
    print ("")

    #Measure divergence from 3D
    
    #--------------------
    #   Lambda methods
    #--------------------
    lambd.lambdacompare(thiscluster) #compare lambda methods (uses 3D)
    filehandling.mergefiles()
    print ("")
    
    #--------------------
    #    k comparisons
    #--------------------
    #compare simulations (k##)
    lambd.lambda_k(thiscluster) #compare simulations (k## together)


#print (timeit.default_timer() - start_time)
