#!/usr/bin/env python

from __future__ import division #import division from V3 to avoid divide by 0
import os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import container #needed to remove error bars from legend
from matplotlib.ticker import AutoMinorLocator #for minor ticks

import argparse

parser = argparse.ArgumentParser(description='Lambda for median MST')

parser.add_argument('model_info',nargs='+')

parser.add_argument('--medians3d', nargs='*', type=float, default=[])
parser.add_argument('--medians2d', nargs='*', type=float, default=[])

parser.add_argument('--quartiles3d', nargs=5, type=float)
parser.add_argument('--quartiles2d', nargs=5, type=float)

parser.add_argument('--sixths3d', nargs=7, type=float)
parser.add_argument('--sixths2d', nargs=7, type=float)

args = parser.parse_args()

plotsdir = args.model_info[0]   # or sys.argv[1]
kval = args.model_info[1]      # or sys.argv[2]
ctype = args.model_info[2].split("_")[1] #remove 'cluster_' prefix

medians3d = args.medians3d
medians2d = args.medians2d
nsnaps = len(args.medians3d)

thissnap = int(args.model_info[3])

#print args.quartiles3d
#print args.quartiles2d

#print plotsdir
#print kval
#print ctype

mpl.rcParams['lines.linewidth'] = 1.0 #set default line width to 1.0
#
projections=['3D','xy']
duration = 10 #duration in Myr (x-axis limit)

tmin = duration/nsnaps
#add extra element to max value as 'arange' goes to max+1:
time = np.arange(tmin, duration + tmin, duration/nsnaps)

#Make plot:
fig, ax = plt.subplots()
#plt.figure()
plt.xlabel("Time (Myr)")
plt.ylabel('$\overline{\Lambda}$')
plt.title("Mass segregation")
#n-1 minor ticks:
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
#plt.xlim(0,duration)
#plt.ylim(0,8)
if time[thissnap-1] < 5:
    xmin = 0 ; xmax = 6
else:
    xmin = 4 ; xmax = 10
plt.xlim(xmin,xmax)
plt.ylim(0,4)

#print thissnap
#print(time[thissnap-1],medians3d[thissnap-1],
#      args.quartiles3d[1],args.quartiles3d[3])


ax.plot(time,medians3d,label='3D')
ax.plot(time,medians2d,label='xy')

#reset colour cycle:
plt.gca().set_prop_cycle(None)

plt.plot(time[thissnap-1],medians3d[thissnap-1],marker='.')
plt.plot(time[thissnap-1],medians2d[thissnap-1],marker='.')

#reset colour cycle:
plt.gca().set_prop_cycle(None)

#error bars at quartiles:
#err3d_lo = medians3d[thissnap-1] - args.quartiles3d[1]
#err3d_hi = args.quartiles3d[3] - medians3d[thissnap-1]
#err2d_lo = medians2d[thissnap-1] - args.quartiles2d[1]
#err2d_hi = args.quartiles2d[3] - medians2d[thissnap-1]
#saveplot = ("%s/snap%04d_lambda_quartiles.png" % (plotsdir,thissnap))
#error bars at 1/6, 5/6:
err3d_lo = medians3d[thissnap-1] - args.sixths3d[1]
err3d_hi = args.sixths3d[-2] - medians3d[thissnap-1]
err2d_lo = medians2d[thissnap-1] - args.sixths2d[1]
err2d_hi = args.sixths2d[-2] - medians2d[thissnap-1]
saveplot = ("%s/snap%04d_lambda_17-83pc.png" % (plotsdir,thissnap))

#shift error bars by tiny amount either side to avoid overlap
errshift = 0.002*(xmax-xmin)

plt.errorbar(time[thissnap-1]+errshift,medians3d[thissnap-1],marker='None',
             yerr=[[err3d_lo],[err3d_hi]],
             capsize=4
)[-1][0].set_linestyle('--')

plt.errorbar(time[thissnap-1]-errshift,medians2d[thissnap-1],marker='None',
             yerr=[[err2d_lo],[err2d_hi]],
             capsize=4
)[-1][0].set_linestyle('--')

plt.legend(loc="upper right", fontsize=10)

plt.tight_layout() #smaller margins
#plt.show()
plt.savefig(saveplot, bbox_inches='tight')
plt.close()
