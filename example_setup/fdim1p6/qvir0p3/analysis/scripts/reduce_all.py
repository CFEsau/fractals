#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import numpy as np
import sys

ic  = str(raw_input('input the path to the directory with the .sl run files: '))
out  = str(raw_input('input the path to the output directory: ')) 
    
nfile=0
for fname in os.listdir(ic + '/'):
    if 'run' in fname and '.sl' in fname:
        nfile += 1
        print nfile
        os.system('../reduce ' + ic + '/' + fname + ' ' + out + '/' + fname[0:-3])

print 'Number of files: ' + str(nfile)
print 'Done'
