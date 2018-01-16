#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import numpy as np
import sys

fbinval=''
while not fbinval:
    fbinval=str(raw_input('Enter binary fraction: '))
    if fbinval == '0%':
        fbin='fbinary0p0'
        break
    elif fbinval == '50%':
        fbin='fbinary0p5'
        break
    elif fbinval == '100%':
        fbin='fbinary1p0'
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
        fdim='f16'
        break
    elif fdimval == '2.0':
        fdim='f20'
        break
    elif fdimval == '2.6':
        fdim='f26'
        break
    elif fdimval == '3.0':
        fdim='f30'
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
        qvir='q03'
        break
    elif qvirval == '0.5':
        qvir='q05'
        break
    else:
	print "\n Warning: virial ratio %s not recognised." % qvirval
	print "     Options are 0.3, 0.5."
	print "     Enter new virial ratio: "
        qvirval=''

#modelpath='../'+fbin+'/'+fdim+qvir
modelpath='/local/cfe/backed_up_on_astro3/fractals/r1p0/'+fbin+'/'+fdim+qvir

ic  = str(raw_input('input the name of the directory with the .sl run files: '))
outdir  = str(raw_input('input the name of the output directory: '))

nfile=0
for fname in os.listdir(modelpath+'/'+ic +'/'):
    if 'run' in fname and '.sl' in fname:
        nfile += 1
        print nfile
        #run 'reduce' executable:
        os.system('./reduce ' +modelpath+'/'+ic+'/'+fname+ ' '
                  + modelpath+'/'+outdir+'/'+fname[0:-3])

print 'Number of files: ' + str(nfile)
print 'Done'
