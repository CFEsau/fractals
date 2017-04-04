#!/usr/bin/env python

from __future__ import print_function
import os,sys,time
import plotconfig,filehandling
import energy,virial,rhalfm,lambd
#import timeit

skip='y'                   #skip energy, Q, rhalf plots? y/n
plotconfig.errorbars='n'   #plot error bars on lambda?
plotconfig.duration = 10   #Duration of simulation (Myr)
    
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
    energy.printenergy()
    energy.plotenergy()
    filehandling.mergefiles()

    virial.printvirial()
    virial.plotvirial()
    filehandling.mergefiles()

    rhalfm.printrhalfm()
    #do for different types of cluster:
    #print ("Cluster list:",plotconfig.clustertypes)
    for thiscluster in plotconfig.clustertypes:
        rhalfm.plotrhalfm(thiscluster)
        filehandling.mergefiles()
        print ("")
    
else:
    print("Skipping energy...")
    print("Skipping virial ratio...")
    print("Skipping half-mass radius...")

    
lambd.printlambda()
lambd.lambdasetup()

#start_time = timeit.default_timer()
#overplot projections for all types of lambda, for each cluster

#compare projections for each lambda
for thiscluster in plotconfig.clustertypes:
    lambd.lambdaprojections(thiscluster)
    filehandling.mergefiles()
    print ("")

#print (timeit.default_timer() - start_time)

#compare lambda methods (uses 3D projection)
for thiscluster in plotconfig.clustertypes:
    lambd.lambdacompare(thiscluster)
    filehandling.mergefiles()
    print ("")

#compare simulations (k##)

'''
print ("")
print ("************************************")
print ("* Making comparison plots across k *")
print ("************************************")
# Plot all values of k together for each E, Q, lambda, etc

if skip == 'n':
    energy.energy_k()
    virial.virial_k()
    for thiscluster in plotconfig.clustertypes:
        rhalfm.rhalfm_k(thiscluster)
        
#for thiscluster in plotconfig.clustertypes:
    #lambd.lambda_k(thiscluster)
'''
'''
#original bash script:


echo ""
echo "****************************"
echo "* Comparing lambda methods *"
echo "****************************"
# Plot lambar, lamtilde, etc against each other for each k

plots/differentlambda.py ${fbin} ${fdim} ${qvir} ${outpath}

#param=alllam
. ./plots/merge.sh

echo ""
echo "       ... done."


echo ""
echo "************************************"
echo "* Making comparison plots across k *"
echo "************************************"
# Plot all values of k together for each E, Q, lambda, etc

plots/lambda_comparek.py ${fbin} ${fdim} ${qvir} ${outpath}
echo ""
sleep 0.6

echo "       ... done."
'''
