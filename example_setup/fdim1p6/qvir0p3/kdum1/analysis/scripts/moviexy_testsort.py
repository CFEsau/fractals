#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
import numpy as np
import sys

my_dpi=96
fig = plt.figure(figsize=(960/my_dpi, 960/my_dpi), dpi=my_dpi)
arglist = sys.argv
sim = ''

if(len(arglist) < 2):
    #sim  = str(raw_input('input the path to the simulation directory: '))
    sim = 'analysis_test/outputs/runinv_0100/'
else:
    sim = str(arglist[1])

print sim

nfile=0
for fname in os.listdir(sim + '/snapshots/'):
    if 'snap' in fname and not 'tmp' in fname and not 'png' in fname:
        nfile += 1
print 'Number of files: ' + str(nfile)

for i in range(1,nfile+1):

    ifname = sim + '/snapshots' + '/snap'
    if i < 10:
        ifname += '000' + str(i)
    elif i < 100:
        ifname += '00' + str(i)
    elif i < 1000:
        ifname += '0' + str(i)
    elif i < 10000:
        ifname += '' + str(i)
    else:
        ifname += str(i)

    print ifname
    
    fo = open(ifname, "r")
    lines = fo.readlines();
    i=0
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
    m=[]
    numbers=[]
    nstars=0
#for finding the n most massive stars in the system:
    most_massive=[]
    n_list=5 #number of stars to go in list

    for line in lines:
        numbers.append(line.split())
        x.append(np.float64(numbers[nstars][4]))
        y.append(np.float64(numbers[nstars][5]))
        z.append(np.float64(numbers[nstars][6]))
        m.append(np.float64(numbers[nstars][3]))

        nstars+=1

    if ifname==sim+'/snapshots'+'/snap0001':
        most_massive=sorted(m)
        most_massive=most_massive[-n_list:]
        print 'snap0001:',most_massive
        print 'final in list: ',most_massive[-1]
    if ifname==sim+'/snapshots'+'/snap0099':
        most_massive=sorted(m)
        most_massive=most_massive[-n_list-1:]
        print 'snap0099:',most_massive

    fo.close()

    plt.clf()
    plt.locator_params(nbins=16)
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y,marker='.',s=2)
    plt.xlabel('X-Position (pc)')
    plt.ylabel('Y-Position (pc)')
    plt.title('Simulation')
    plt.grid(True)
    plt.axis([-5, 5, -5, 5])
    plt.savefig(ifname+'.xy.png', dpi=my_dpi, bbox_inches='tight')

os.system('avconv -y -r 15 -i ' + sim + '/snapshots/snap%04d.xy.png -s 1024x800 ' + sim + '/xy_test.mp4')

for fn in os.listdir(sim + '/snapshots/'):
    if '.xy.png' in fn:
        os.system("rm -rf " + sim + '/snapshots/' + fn)
