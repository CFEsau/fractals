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
    nth_massive=100000. #nth most massive star - least massive in list
    n_massive=5 #number of stars to go in list

    for line in lines:
        numbers.append(line.split())
        x.append(np.float64(numbers[nstars][4]))
        y.append(np.float64(numbers[nstars][5]))
        z.append(np.float64(numbers[nstars][6]))
        m.append(np.float64(numbers[nstars][3]))

        if ifname==sim+'/snapshots'+'/snap0001':
#============================#
#Finding n most massive stars#
#============================#

            #mass of star at current position in list
            this_m= float(numbers[nstars][3])
        
            #nstars is less than number required in list:
            if nstars<n_massive:                
                #add this_m to the list
                most_massive.append(this_m)
               # print 'this_m:', this_m, 'nth massive:', nth_massive
                #if this_m less than least massive in list, save as nth_massive
                if this_m<nth_massive:
                    nth_massive=this_m
                    
            #nstars is greater than number requried in list:
            else:
                if this_m>nth_massive:
                    #overwrite smallest value in list (at 0 from sorting)
                    most_massive[0]=this_m
                    most_massive=sorted(most_massive)
                    #set nth_massive as first element in sorted list:
                    nth_massive=most_massive[0]

#============================#

        nstars+=1

    if ifname==sim+'/snapshots'+'/snap0001': 
        print 'After sorting: ',most_massive

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
