from __future__ import print_function
import glob,os,time
from PyPDF2 import PdfFileMerger
import plotconfig

def getclusters():
    #directory format: outpath/runinv_k/cluster/.dat
    #can't use 1st k as there might be a later k with different cluster.
    #use 1st to make list & look through others to add to list?
    #then when looping through, have an 'if doesn't exist' statement
    print("")
    owd = os.getcwd()
    os.chdir(plotconfig.outpath)
    for kmodel in os.listdir('./'):
        #print ("This model:",kmodel)
        if 'runinv_k' in kmodel:
            for thiscluster in os.listdir(kmodel):
                if 'cluster_' in thiscluster:
                    if thiscluster not in plotconfig.clustertypes:
                        plotconfig.clustertypes.append(thiscluster)
    os.chdir(owd)
    

def mergefiles():
    #
    # Go into 'plots' directory for this fdim, qvir.
    # List files beginning with 'k##_...'
    # List all files with pattern 'k##_...parameter...'
    # Merge these files & save under 'parameter_cluster.pdf'
    #file format: {filepath}/k#_{cluster}_{param}.pdf
    #
    owd = os.getcwd()
    filepath = plotconfig.outpath+'/plots'
    os.chdir(filepath)

    flist = glob.glob('k*.pdf') 
       
    while flist:  #do while flist isn't empty
        thiskplot=flist[0]  #take first filename
        filestructure = thiskplot.split('_')
        #kval = filestructure[0]
        thiscluster = filestructure[1]
        parameter = (filestructure[2]).split('.',1)[0] #strip .pdf extension
        
        pdflist = sorted(glob.glob('k*_'+thiscluster+'_'+parameter+'.pdf'))
        if pdflist:
            mergelist = PdfFileMerger()
            for thispdf in pdflist:
                mergelist.append(thispdf)
            newfile = parameter+'_'+thiscluster+'.pdf'
            mergelist.write(newfile)
            print ("Files merged: "+newfile)
            for thispdf in pdflist:
                os.remove(thispdf)
            pdflist=[]
        else:
            break

        flist = glob.glob('k*.pdf') #update file list
            
    os.chdir(owd)
    #time.sleep(0.4)
