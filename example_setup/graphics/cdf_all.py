#!/usr/bin/env python

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
#from sys import argv
#outpath=argv[1] # path to 'outputs' directory, passed from cdfplots.sh
from matplotlib.ticker import AutoMinorLocator #for minor ticks
import glob #import 'global module' for computer access
import subprocess #for calling bash script
import itertools #for cycling through different markers in plots (used in Delta)


filename = 'MSTedgeL_3D' #do for 3D simulation (don't include suffix)
cluster='cluster_FoV5pc'
nstars=10 #number of stars in MST
nedge=nstars-1 #number of edge connections in MST

nmst=20 #number of random MSTs

simnum='01' #k-value (1-10, from models)
snapshots=[76, 105, 135, 459] #snapshot number

fbinval=''
while not fbinval:
    fbinval=str(raw_input('Enter binary fraction: '))
    if fbinval == '0%':
        fbin='fbin0p0'
        break
    elif fbinval == '50%':
        fbin='fbin0p5'
        break
    elif fbinval == '100%':
        fbin='fbin1p0'
        break
    else:
        print '\n Warning: Binary fraction %s not recognised.' % fbinval
	print '     Options are 0%, 50%, 100%.'
	print '     Enter new binary fraction: '
        fbinval=''
    
fdimval=''
while not fdimval:
    fdimval=str(raw_input('Enter fractal dimension: '))
    if fdimval == '1.6':
        fdim='fdim1p6'
        break
    elif fdimval == '2.0':
        fdim='fdim2p0'
        break
    elif fdimval == '2.6':
        fdim='fdim2p6'
        break
    elif fdimval == '3.0':
        fdim='fdim3p0'
        break
    else:
        print '\n Warning: fractal dimension %s not recognised.' % fdimval
	print '     Options are 1.6, 2.0, 2.6, 3.0.'
	print '     Enter new fractal dimension: '
        fdimval=''
        
qvirval=''
while not qvirval:
    qvirval  = str(raw_input('Enter virial ratio: '))
    if qvirval == '0.3':
        qvir='qvir0p3'
        break
    elif qvirval == '0.5':
        qvir='qvir0p5'
        break
    else:
	print "\n Warning: virial ratio %s not recognised." % qvirval
	print "     Options are 0.3, 0.5."
	print "     Enter new virial ratio: "
        qvirval=''


outpath='../'+fbin+'/'+fdim+'/'+qvir+'/outputs'
simpath=outpath+'/runinv_k'+simnum+'/'+cluster+'/' #path to data
plotsdir=outpath+'/plots/cdf_k'+simnum #directory for saving plots

#make directory for CDFs for this k (if it doesn't exist)
os.system('mkdir -p '+plotsdir)

#data:
objmacro = np.loadtxt(simpath+'CDFdata/'+filename+'.dat')
#print objmacro
#set minor tick marks:
minorLocator = AutoMinorLocator(4) #this gives n-1 minor ticks

#######################################################################

#plot the cdf for this snapshot
#get probabilities from 0 to 1 for y axis:
probrange=np.linspace(0,1,num=nstars)

#for Delta 1 vs Delta 2 axes:
axislimdelta=0.
marker = itertools.cycle(('+','*', 'x', '.')) 

for thissnap in snapshots:
    print ""
    print "*************"
    print "Snapshot %d:" % thissnap
    print "*************"
    print "----------------"
    print "-- Object CDF --"
    print "----------------"

    objedgelengths=[]
    objedgelengths=objmacro[thissnap-1][1:]
    
    #---------------------------------------
    # Average edge lengths for object stars
    #---------------------------------------
    
    # Edge length data read in above as 'objmacro'
    mean=0.
    medianedge=0.
    median=0.
    geometriclength=0.
    geometric=0.
    
    # Find mean
    treelength=0.
    treelength=sum(objedgelengths)
    mean=treelength/float(nedge)
    
    # Find median
    medianedge=(nedge/2.0)+0.5
    if ((medianedge).is_integer()):
        medianedge=int(medianedge)
        median=objedgelengths[medianedge-1] #-1 as python starts at 0
    else:
        #Take mean of the two median values
        medianedge=int(medianedge-0.5)
        median=objedgelengths[medianedge-1]+objedgelengths[medianedge]
        median=median*0.5
        
    # Find geometric mean
    #Sum the logs of individual lengths
    geometriclength=sum(np.log(objedgelengths))
    geometric=np.exp(geometriclength/float(nedge))
    
    #--------------------------
    # Get x & y values to plot
    #--------------------------
    
    #x-values:
    #Each edge length is recorded twice as the CDF here is a step function
    xedgelengths=[]
    for l in range(0,nedge): #from 1 to nedge in python
        xedgelengths.append(objedgelengths[l])
        xedgelengths.append(objedgelengths[l])

    #set max x-axis limit a little greater than max edge length:
    xmax = 1.05*xedgelengths[2*(nedge-1)]
    #and use this at end of edge lengths list to get horizontal line:
    xedgelengths.append(xmax)
    
    
    #y-values:    
    #Again, CDF is a step function so need these values twice:
    probability=[]
    #1st y value is 0
    probability.append(0)
    #following y values are in the range 0>y>=1:
    for p in range(1,nedge+1): #start from 1 as already got 0th entry
        probability.append(probrange[p])
        probability.append(probrange[p])
    
    #--------------------------------------------------------
    # Get coordinates for markers showing different averages
    #--------------------------------------------------------
    
    # draw vertical lines from (avtype,bottom) to (avtype,top)
    # for each type of average. Iterate through edge lengths
    # until we find where each marker needs to be placed.

    domean=True
    domedian=True
    dogeometric=True
    i=0
    
    while domean:
        if mean<xedgelengths[i]: #find when arithmetic mean < edge L
            topmean=probability[i]+0.05
            botmean=probability[i]-0.05
            xlabelmean=mean-(0.03*xmax) #labels recalculated later for...
            ylabelmean=probability[i]+0.03 #...plots with random stars
            i=0
            domean = False
        else:
            i=i+2 #+2 as values are repeated in list for step function
            
    while domedian:
        if median<xedgelengths[i]: #find when median < edge L
            topmed=probability[i]+0.05
            botmed=probability[i]-0.05
            xlabelmed=median-(0.03*xmax)
            ylabelmed=probability[i]+0.03
            i=0
            domedian = False
        else:
            i=i+2 #+2 as values are repeated in list for step function
            
    while dogeometric:
        if geometric<xedgelengths[i]: #find when geometric mean < edge L
            topgeo=probability[i]+0.05
            botgeo=probability[i]-0.05
            xlabelgeo=geometric-(0.03*xmax)
            ylabelgeo=probability[i]+0.03
            i=0
            dogeometric = False
        else:
            i=i+2 #+2 as values are repeated in list for step function

    
    #-----------------
    # Plot object CDF
    #-----------------

    #Do one plot for just the 'object' CDF:
    my_dpi=96
    saveplot = ''
    plt.figure()
    plt.xlabel("Edge length (pc)")
    plt.ylabel("Probability")
    plt.title("Edge lengths of object stars")
    plt.xlim(-0.01,xmax)
    plt.ylim(ymax = 1.005, ymin = 0)
    plt.axes().xaxis.set_minor_locator(minorLocator)
    
    #plot data:    
    plt.plot(xedgelengths,probability)
    
    plt.annotate("snap "+str(thissnap), xy=(0.95,-0.09),
                 xycoords='axes fraction', fontsize = 11)
    
    #plot markers showing different types of mean for object stars:
    plt.plot([mean,mean], [botmean,topmean], 'k--', lw=0.7)
    plt.plot([median,median], [botmed,topmed], 'r--', lw=0.7)
    plt.plot([geometric,geometric], [botgeo,topgeo], 'g--', lw=0.7)
    #and labels:
    plt.text(xlabelmean, ylabelmean, '$\overline{\Lambda}$',color='k',fontsize=11)
    plt.text(xlabelmed, ylabelmed, '$\widetilde{\Lambda}$',color='r',fontsize=11)
    plt.text(xlabelgeo, ylabelgeo, '$\Gamma$',color='g',fontsize=11)
    
    #plt.show()
    saveplot = (plotsdir+'/cdf_snap'+'%03d'%thissnap+'_obj.png')
    plt.savefig(saveplot, bbox_inches='tight')
    print "Plot saved at " + saveplot
    plt.close()
    
    
#End of 'object' CDF
    
    #######################################################################
    print ""
    print "-----------------"
    print "-- Random CDFs --"
    print "-----------------"
    
    #allranedgelengths=np.zeros(shape=(nedge,nmst),order='F')
    
    #create array of nedge by nmst: (list of lists)
    #(could do this with one list & put in same loop below but need all
    #data if you want to have all x axes on the same scale across CDFs)

    xmaxran=0. #for maximum x value calculated in loop
    allranedgelengths=[]
    
    for thisMST in range(1,nmst+1): #i=1,2,...20 for nmst=20
        ranedgelengths=[]
        
        #load edge lengths file (edge lengths for all snapshots)
        ranmacro=np.loadtxt(simpath+'CDFdata/'+filename+'_'+str(thisMST)+'.dat')
        #list of edge lengths for given snapshot without snapshot# at beginning
        ranedgelengths=ranmacro[thissnap-1,1:]

        #Find maximum edge length across all CDFs:
        if ranedgelengths[-1]>xmaxran:
            xmaxran=ranedgelengths[-1]
            
        #make nested list of all edge lengths for each MST:
        #(needed for Delta1 vs Delta2 plot)
        allranedgelengths.append(ranedgelengths[:])
        #[[MST1L1,MST1L2,MST1L3...],[MST2L1,MST2L2,MST2L3...],...]
        #referenced as allranedgelengths[mst#-1][edgelength#-1]
        
        
    #---------------------------------------
    # Average edge lengths for random stars
    #---------------------------------------
    #Find average edge lengths & label positions and plot each CDF:
    for thisMST in range(1,nmst+1):
        #must be a separate loop as we needed to run
        #through them all before to find xmax.
        
        ranedgelengths=[]
        ranedgelengths=allranedgelengths[thisMST-1][:]
        
        # Edge length data for all CDFs is in ranedgelengths[nmst][nedge]
        # (remember lengths are saved twice for step function)
        
        treelength=0.
        meanran=0.
        medianedge=0.
        medianran=0.
        geometriclength=0.
        geomeanran=0.
        
        # Find mean
        #sum edges for total tree length:
        treelength=sum(ranedgelengths)
        #treelength=sum(objmacro[thissnap-1,1:nstars])
        meanran=treelength/float(nedge)
        
        # Find median
        medianedge=(nedge*0.5)+0.5
        if ((medianedge).is_integer()):
            medianedge=int(medianedge)
            medianran=ranedgelengths[medianedge-1] #-1 as python starts at 0
        else:
            #Take mean of the two median values
            medianedge=int(medianedge-0.5)
            medianran=ranedgelengths[medianedge-1]+ranedgelengths[medianedge]
            medianran=medianran*0.5
           
        # Find geometric mean
        #Sum the logs of individual lengths
        geometriclength=sum(np.log(ranedgelengths))
        geomeanran=np.exp(geometriclength/float(nedge))
    
        #----------------------
        # Get x values to plot
        #----------------------
        # y same as 'object' values
        
        #Each edge length is recorded twice as the CDF is a step function:
        xlengthsran=[]
        for l in range(0,nedge):
            xlengthsran.append(ranedgelengths[l])
            xlengthsran.append(ranedgelengths[l])        
        #add slightly larger value than xmaxran for last entry:
        xlengthsran.append(1.05*xmaxran)
        
        #--------------------------------------------------------
        # Get coordinates for markers showing different averages
        #--------------------------------------------------------
        
        # draw vertical lines from (avtyperan,bottom) to (avtyperan,top)
        # for each type of average. Iterate through edge lengths
        # until we find where each marker needs to be placed.
        
        domean=True
        domedian=True
        dogeometric=True
        i=0
        
        while domean:
            #find first instance of mean being less:
            if meanran<xlengthsran[i]:
                topmeanran=probability[i]+0.05
                botmeanran=probability[i]-0.05
                xlabelmeanran=meanran-(0.03*xmax)
                ylabelmeanran=probability[i]+0.03
                i=0
                domean = False
            else:
                i=i+2 #+2, values repeated in list for step function
                
        while domedian:
            #find first instance of median being less
            if medianran<xlengthsran[i]:
                topmedran=probability[i]+0.05
                botmedran=probability[i]-0.05
                xlabelmedran=medianran-(0.03*xmax)
                ylabelmedran=probability[i]+0.03
                i=0
                domedian = False
            else:
                i=i+2 #+2, values repeated in list for step function
            
        while dogeometric:
            #find first instance of geometric mean being less
            if geomeanran<xlengthsran[i]:
                topgeoran=probability[i]+0.05
                botgeoran=probability[i]-0.05
                xlabelgeoran=geomeanran-(0.03*xmax)
                ylabelgeoran=probability[i]+0.03
                i=0
                dogeometric = False
            else:
                i=i+2 #+2, values repeated in list for step function

                
        #---------------------------
        # Plot CDFs for random MSTs
        #---------------------------
        
        my_dpi=96
        saveplot=''
        plt.figure()
        plt.xlabel("Edge length (pc)")
        plt.ylabel("Probability")
        plt.title("Edge lengths of random MSTs")
        xmaxindi=max(ranedgelengths[-1],objedgelengths[-1])
        plt.xlim(-0.01,1.05*xmaxindi) #for individual xmax
        #plt.xlim(-0.01,1.05*xmaxran) #for global xmax
        plt.ylim(ymax = 1.005, ymin = 0)
        plt.axes().xaxis.set_minor_locator(minorLocator)
        #plt.axes().yaxis.set_minor_locator(minorLocator*0.5)

        #plot data:
        plt.plot(xlengthsran,probability)
        #plot markers showing different types of mean for object stars:
        plt.plot([meanran,meanran], [botmeanran,topmeanran], 'k--', lw=0.7)
        plt.plot([medianran,medianran], [botmedran,topmedran], 'r--', lw=0.7)
        plt.plot([geomeanran,geomeanran], [botgeoran,topgeoran], 'g--', lw=0.7)
        #and labels:
        plt.text(xlabelmeanran, ylabelmeanran, '$\overline{\Lambda}$',color='k',fontsize=11)
        plt.text(xlabelmedran, ylabelmedran, '$\widetilde{\Lambda}$',color='r',fontsize=11)
        plt.text(xlabelgeoran, ylabelgeoran, '$\Gamma$',color='g',fontsize=11)
        #---------------------
        #End of 'random' CDFs
        #---------------------
        
        #plot 'object' cdfs for comparison:
        plt.plot(xedgelengths,probability,color='grey')
        #plot markers showing different types of mean for object stars:
        plt.plot([mean,mean], [botmean,topmean], 'k--', lw=0.7)
        plt.plot([median,median], [botmed,topmed], 'r--', lw=0.7)
        plt.plot([geometric,geometric], [botgeo,topgeo], 'g--', lw=0.7)
        #and labels:
        #plt.text(xlabelmean, ylabelmean, '$\overline{\Lambda}$',color='k',fontsize=11,alpha=0.5)
        #plt.text(xlabelmed, ylabelmed, '$\widetilde{\Lambda}$',color='r',fontsize=11,alpha=0.5)
        #plt.text(xlabelgeo, ylabelgeo, '$\Gamma$',color='g',fontsize=11,alpha=0.5)
        #----------------------------
        #End of greyed 'object' CDFs
        #----------------------------

        
        plt.annotate("snap "+str(thissnap), xy=(0.86,0.03),
                     xycoords='axes fraction', fontsize = 11)
        #plt.show()
        saveplot = (plotsdir+'/cdf_snap'+'%03d'%thissnap+'_'+'%02d'%thisMST+'.png')
        plt.savefig(saveplot, bbox_inches='tight')
        print "Plot saved at " + saveplot
        plt.close()

        
    #Merge CDFs for this snapshot into one document:
    print ""
    print "     Merging CDF plots..."
    filestring=str("")
    for files in glob.glob(plotsdir+'/cdf_snap'+'%03d'%thissnap+'_*.png'):
        filestring=filestring+str(" ")+files
        #print filestring
    os.system('convert '+filestring+' '+plotsdir+'/cdf_snap'+'%03d'%thissnap+'.pdf')
    print ""
    print "                Saved in cdf_snap"+'%03d'%thissnap+".pdf"
    print ""
    #remove snapshots:
    os.system('rm '+filestring)
    #End of 'random' CDFs
    
    
    #######################################################################
    print ""
    print "--------------"
    print "-- D1 vs D2 --"
    print "--------------"
    #for each snapshot plot scatter of Delta1 vs Delta2 for each CDF
    #Delta is the difference between the object star edge length and
    #random MST edge length, so there will be nmst values per plot.
    
    edgenum=[3,7]

    #Get object star edge lengths at edgenum[0] & edgenum[1]:
    deltaobj=[0,0]
    deltaobj[0]=objedgelengths[edgenum[0]-1]#want 3rd in list, so py index=2
    deltaobj[1]=objedgelengths[edgenum[1]-1]
    
    #Get random star edge lengths at edgenum[0] & edgenum[1]:
    deltax=[] #nmst entries
    deltay=[]
    #make empty array 'deltapoints' size 2xnmst
    deltapoints=[]
    for i in range (2):
        deltapoints.append([])
        for j in range (0,nmst):
            deltapoints[i].append(0)
    
    for thisMST in range (1,nmst+1):
        deltaran=[0,0]
        deltaran[0]=allranedgelengths[thisMST-1][edgenum[0]-1]
        deltaran[1]=allranedgelengths[thisMST-1][edgenum[1]-1]
        
        #take difference between edge lengths:
        deltapoints[0][thisMST-1]=deltaran[0]-deltaobj[0]
        deltapoints[1][thisMST-1]=deltaran[1]-deltaobj[1]
        
    # -----------------
    # plot dx1 vs dx2
    #------------------

    # all snapshots on same plot:
    saveplot0=''
    plt.figure(0)
    plt.xlabel("$\Delta$"+str(edgenum[0])+"$_{\mathrm{ran-obj}}$")
    plt.ylabel("$\Delta$"+str(edgenum[1])+"$_{\mathrm{ran-obj}}$")
    #find max dx value for axis limits:
    axislimdelta = max(np.max(deltapoints),np.min(deltapoints),
                       axislimdelta,key=abs)
    #if xmax or |xmin| is greater than before, the limits will update
    plt.xlim(-axislimdelta*1.05,axislimdelta*1.05)
    plt.ylim(-axislimdelta*1.05,axislimdelta*1.05)
    
    #plot data:
    print "    Populating delta plot for snapshot " + str(thissnap)
    plt.plot(deltapoints[0][:],deltapoints[1][:],linestyle="",marker=marker.next(),
             label=thissnap)
    
    plt.legend(fontsize=10,title='snapshot #',numpoints=1)
    #plt.show(0)
    #----------------------------

    #snapshots on separate plots:
    #savethissnap=''
    plt.figure(thissnap)
    #plt.xlabel("$\Delta$1")
    #plt.ylabel("$\Delta$2")
    #find max dx value for axis limits:
    #axislim = max(np.max(deltapoints),np.min(deltapoints),
                  #axislimdelta,key=abs)
    #plt.xlim(-axislim*1.05,axislim*1.05)
    #plt.ylim(-axislim*1.05,axislim*1.05)
    #plt.annotate("snaphot " + str(thissnap), xy=(-0.06,-0.09),
    #             xycoords='axes fraction', fontsize = 10)
    
    #plot data:
    plt.plot(deltapoints[0][:],deltapoints[1][:],linestyle="",marker="x")
    
    #plt.show()
    #savethissnap.append(plotsdir+'/snap'+'%03d'%thissnap+'_delta.png')
    #plt.savefig(saveplotthissnap, bbox_inches='tight')
    #print "Graph saved at " + savethissnap
    #plt.close(thissnap)

#
#End of 'snapshot' loop

#save & close plot 0:
plt.figure(0)
saveplot0 = (plotsdir+'/snapall_delta.png')
plt.savefig(saveplot0, bbox_inches='tight')
print "Plot saved at " + saveplot0

#Do individual D1 by D2 plots (need separate loop for same axis lims):
for thissnap in snapshots:
    savesnap=''
    plt.figure(thissnap)
    plt.xlabel("$\Delta$"+str(edgenum[0])+"$_{\mathrm{ran-obj}}$")
    plt.ylabel("$\Delta$"+str(edgenum[1])+"$_{\mathrm{ran-obj}}$")
    plt.xlim(-axislimdelta*1.05,axislimdelta*1.05)
    plt.ylim(-axislimdelta*1.05,axislimdelta*1.05)
    plt.annotate("snaphot " + str(thissnap), xy=(-0.06,-0.09),
                 xycoords='axes fraction', fontsize = 10)
    #plt.plot(deltapoints[0][:],deltapoints[1][:],linestyle="",marker="x")
    savesnap=(plotsdir+'/snap'+'%03d'%thissnap+'_delta.png')
    plt.savefig(savesnap, bbox_inches='tight')
    print "D1 vs D2 saved at " + savesnap


#Merge D1/D2 for each snapshot into one document:
print ""
print "   Merging Delta plots..."
filestring=str("")
for files in glob.glob(plotsdir+'/snap*_delta.png'):
    filestring=filestring+str(" ")+files
os.system('convert '+filestring+' '+plotsdir+'/deltaedges.pdf')
print ""
print "           deltaedges.pdf"
print ""
#remove plots for individual snapshots:
os.system('rm '+filestring)
#

#Plot all D1/D2 on one plot:
#Need two nested lists, one for D1 & one for D2.
#[[x1,x2...x10]snap1,[x1...x10]snap2...,[x...]snapmax],[[same for y],[],...[]]
#OR find way to have two plots on the go at once so global plot can be
#made in same loop as other D1/D2 plot, and set axis lims at the end.
#######################################################################

#---------------------------
# Plot lambdas with snap #s
#---------------------------

print ""
print "------------------"
print "-- Lambda plots --"
print "------------------"

duration = 10. #duration of simulation (Myr)

plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel(r"$\Lambda$")
plt.title("Mass segregation")
plt.ylim(0,25)
plt.axes().xaxis.set_minor_locator(minorLocator)

#plot different lambda in certain order:
#(not bothering with y/n like befoer...
#to add more look at differentlambda.py)
lambda_types=[]
lambda_types.append('lambar')
lambda_types.append('lamrms')
lambda_types.append('lamtil')
lambda_types.append('gam')
#dictionary for lambda data:
lambda_data=dict()

#dictionary for tex strings:
lambda_tex=dict()
lambda_tex['lambar']='$\overline{\Lambda}$'
lambda_tex['lamrms']='$\Lambda_{rms}$'
lambda_tex['lamtil']='$\widetilde{\Lambda}$'
lambda_tex['gam']='$\Gamma$'


#load lambda data:
for lamfile in os.listdir(simpath+'lambda/'):
    if '_3D' in lamfile:
        lamtype = lamfile.split("_")[1] #get lambda type (bar, rms, etc)
        lam = np.loadtxt(simpath+'lambda/'+lamfile) #load data
        nsnap =  lam[:,0]
        time = (nsnap/nsnap[-1])*duration
        lambda_data[lamtype] = lam[:,3] # save data for each lamtype
        
#plot each lambda type:
for lamtype in lambda_types:
    plt.plot(time,lambda_data[lamtype],label=lambda_tex[lamtype])
    plt.legend(fontsize=10)
    plt.annotate("fbin = " + fbinval, xy=(0.95,-0.09), 
                 xycoords='axes fraction', fontsize = 10)
    plt.annotate('fdim = ' + fdimval + ', qvir = ' + qvirval, 
                 xy=(-0.06,-0.09), xycoords='axes fraction', fontsize = 10)
    
#draw line at chosen snapshots:
for eachsnap in snapshots:
    thistime=(eachsnap/nsnap[-1])*duration
    plt.axvline(x=thistime,color='k',linestyle='dashed', linewidth=0.5)
        
#plt.show()
saveplot = (plotsdir + '/lambdasnapshots.pdf')
plt.tight_layout()
plt.savefig(saveplot, bbox_inches='tight')
print "   Plot saved at %s" % saveplot
print ""
plt.close()
