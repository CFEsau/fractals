#!/usr/bin/env python

#from __future__ import print_function
import os,sys,time
import plotconfig,filehandling
import energy,virial,rhalfm,lambd
#import timeit

skip = 'y'                 #skip energy, Q, rhalf plots? y/n
plotconfig.errorbars = 'n' #plot error bars on lambda?
plotconfig.duration = 10   #Duration of simulation (Myr - for axis limits)

#Haven't made 'projection' a variable. Change manually. 3D at the moment.

plotconfig.fbin = ''
while not plotconfig.fbin:
    plotconfig.fbin = raw_input('Enter binary fraction: ')
    if plotconfig.fbin == '0%':
        fbinstr = 'fbinary0p0'
    elif plotconfig.fbin == '50%':
        fbinstr = 'fbinary0p5'
    elif plotconfig.fbin == '100%':
        fbinstr = 'fbinary1p0'
    else:
        print "\n Warning: Binary fraction %s not recognised." % (
            plotconfig.fbin)
        print "\tOptions are 0%, 50%, 100%."
        plotconfig.fbin = ''

plotconfig.fdim = ''
while not plotconfig.fdim:
    plotconfig.fdim = raw_input('Enter fractal dimension: ')
    if plotconfig.fdim == '1.6':
        fdimstr = 'f16'
    elif plotconfig.fdim == '2.0':
        fdimstr = 'f20'
    elif plotconfig.fdim == '2.6':
        fdimstr = 'f26'
    elif plotconfig.fdim == '3.0':
        fdimstr = 'f30'
    else:
        print "\n Warning: Fractal dimension %s not recognised." % (
            plotconfig.fdim)
        print "\tOptions are 1.6, 2.0, 2.6, 3.0."
        plotconfig.fdim = ''

plotconfig.qvir = ''
while not plotconfig.qvir:
    plotconfig.qvir = raw_input('Enter virial ratio: ')
    if plotconfig.qvir == '0.3':
        qvirstr = 'q03'
    elif plotconfig.qvir == '0.5':
        qvirstr = 'q05'
    else:
        print "\n Warning: Virial ratio %s not recognised." % (
            plotconfig.qvir)
        print "\tOptions are 0.3, 0.5."
        plotconfig.qvir = ''

#modelpath = "../../%s/%s%s" % (fbinstr, fdimstr, qvirstr)
modelpath = "/local/cfe/backed_up_on_astro3/fractals/r1p0/%s/%s%s"% (
    fbinstr, fdimstr, qvirstr)
print "\nModel path is %s" % modelpath

print "\n**************************"
print "*      Making plots      *"
print "**************************"

plotconfig.outpath = '%s/analysis' % modelpath
print (plotconfig.outpath)

#get different types of cluster (all, FoV, etc):
plotconfig.clustertypes = []
#filehandling.getclusters() # get list of cluster types
plotconfig.clustertypes.append('cluster_FoV5pc')
#print ("Cluster list:",plotconfig.clustertypes)

#Make 'plots' directory if it doesn't exist:
if not os.path.exists(plotconfig.outpath+'/plots'):
    print "Making new directory %s/plots...\n" % plotconfig.outpath
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
        print "\n"
    
else:
    print "Skipping energy..."
    print "Skipping virial ratio..."
    print "Skipping half-mass radius..."

    
lambd.printlambda() #print banner
lambd.lambdasetup() #set which lambda types to plot & y-label text

#start_time = timeit.default_timer()
#overplot projections for all types of lambda, for each cluster

for thiscluster in plotconfig.clustertypes:
    #--------------------
    # Lambda projections
    #--------------------
    #saved as plots/k##_cluster_lambdatype.pdf before merging
    
    lambd.lambdaprojections(thiscluster) #compare projections for each lambda
    filehandling.mergefiles() #merge plots into one document

    #plot 2D projections relative to 3D for each lambda:
    #saved as plots/k##_cluster_lambdatype_proj.pdf before merging
    lambd.projectioncompare(thiscluster)
    filehandling.mergefiles() #merge plots into one document

    #Measure divergence from 3D

    #--------------------
    #   Lambda methods
    #--------------------
    lambd.lambdacompare(thiscluster) #compare lambda methods (uses 3D)
    filehandling.mergefiles()
    print "\n"
    
    #--------------------
    #    k comparisons
    #--------------------
    #compare simulations (k##)
    lambd.lambda_k(thiscluster) #compare simulations (k## together)


#print (timeit.default_timer() - start_time)

print "\nDone",modelpath
