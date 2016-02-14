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
    sim  = str(raw_input('input the path to the simulation directory: '))
else:
    sim = str(arglist[1])

print sim

nfile=0
for fname in os.listdir(sim + '/snapshots/'):
    if 'snap' in fname and not 'tmp' in fname and not 'png' in fname:
        nfile += 1
print 'Number of files: ' + str(nfile)

for i in range(1,nfile+1):

    ifname = sim + '/snapshots/' + 'snap'
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
    cluster_2D=[]
    numbers=[]
    nstars=0
    n_list=int(10)	# This defines the nth largest sample
    time=0. #time of snapshot in Myr
    
    for line in lines:
        numbers.append(line.split()) #This is tuple of whole table
        x = np.append(x,float(numbers[nstars][4]))
        y = np.append(y,float(numbers[nstars][5]))
        z = np.append(z,float(numbers[nstars][6]))
        m = np.append(m,float(numbers[nstars][3]))
        cluster_2D = np.append(cluster_2D, str(numbers[nstars][10]))
        nstars+=1
    
    x_cluster = x[np.where(cluster_2D=="T")]
    y_cluster = y[np.where(cluster_2D=="T")]
    m_cluster = m[np.where(cluster_2D=="T")]

    time=numbers[0][2] #time of snapshot is 1st row (any will do), 2nd col
    #print time
    
    #sort m in ascending order:

    mass_srt_by_m = sorted(m_cluster)
    #select the n most massive stars:
    mass_selected = mass_srt_by_m[-n_list:]
    
    #print the n most massive stars for each snapshot:
    #print mass_selected
    
    #print the ID numbers for each mass
    #print "SORTED", np.where(m>=mass_selected[0])

    #print results. only need to do this once but do twice to be sure:
    if ifname==sim+'/snapshots/'+'snap0001' or ifname==sim+'/snapshots/'+'snap0999':
        print "Time: %.2f Myr" % (round(float(time),2))
        print "   Masses and indices: "
        for index in range(0,n_list):
            #strip unwanted chars from np.where output:
            mass_index=str(np.where(m_cluster==mass_selected[index]))
            mass_index=mass_index.strip('(array[]),')
            #call i in the output file i_index so it runs from 1 up,
            #not 0 up as in mass_index
            i_index=int(mass_index)+1
            print "   ", mass_selected[index], "   ", i_index

#    print x[np.where(m>=mass_selected[0]) and  np.where(cluster_2D=="T")]
    fo.close()

    plt.clf()
    plt.locator_params(nbins=16)
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y, marker='.',color='black', s=2)

    ax1.scatter(x_cluster[np.where(m_cluster>=mass_selected[0])], y_cluster[np.where(m_cluster>=mass_selected[0])], marker='*',color='red', s=20)

    plt.xlabel('X-Position (pc)')
    plt.ylabel('Y-Position (pc)')
    plt.title('Simulation')

    plt.grid(True)
    plt.axis([-5, 5, -5, 5])
    textstring='Time: %.2f Myr' % (round(float(time),2))
    plt.text(3, 4, textstring)
    plt.savefig(ifname+'.xy.png', dpi=my_dpi, bbox_inches='tight')
    del(mass_srt_by_m)
    del(mass_selected)
os.system('avconv -y -r 24 -i ' + sim + '/snapshots/snap%04d.xy.png -s 1024x800 ' + sim + '/xy.mp4')

for fn in os.listdir(sim + '/snapshots/'):
    if '.xy.png' in fn:
        os.system("rm -rf " + sim + '/snapshots/' + fn)
