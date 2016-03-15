#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import pp
import time
#Global Variables#####

picture_count = int(0)
ncores = int(0)
projection=str('2D')

def generate_snapshot(ifname,snapname,projection,qvir,fdim,kval):
        import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
	import numpy as np

	xy_box = int(5)
	z_box = (1)
	my_dpi=96

	fig = plt.figure(figsize=(960/my_dpi, 960/my_dpi), dpi=my_dpi)
		
	fo = open(ifname, "r")
	
	lines = fo.readlines();

	x=[]
	y=[]
	z=[]
	vx=[]
	vy=[]
	vz=[]
	m=[]
	numbers=[]
	cluster_2D =[]
	cluster_3D =[]
	nstars=0
	n_list=int(10)
	# This defines the nth largest sample
	time=0. #time of snapshot in Myr
	for line in lines:
		numbers.append(line.split()) #This is tuple of whole table
		x = np.append(x,float(numbers[nstars][4]))
		y = np.append(y,float(numbers[nstars][5]))
		z = np.append(z,float(numbers[nstars][6]))
		m = np.append(m,float(numbers[nstars][3]))
		cluster_2D = np.append(cluster_2D, str(numbers[nstars][10]))		
		cluster_3D = np.append(cluster_3D, str(numbers[nstars][11]))				
		nstars+=1
	
	if projection == "3D":		
		x_cluster = x[np.where(cluster_3D=="T")]
		y_cluster = y[np.where(cluster_3D=="T")]
		z_cluster = z[np.where(cluster_3D=="T")]		
		m_cluster = m[np.where(cluster_3D=="T")]
	
	elif projection == "2D":		
		x_cluster = x[np.where(cluster_2D=="T")]
		y_cluster = y[np.where(cluster_2D=="T")]
		m_cluster = m[np.where(cluster_2D=="T")]		
	
	time=numbers[0][2] #time of snapshot is 1st row (any will do), 2nd col
	#print time
	mass_srt_by_m = sorted(m_cluster)
	mass_selected = mass_srt_by_m[-n_list:]

	#print "SORTED", np.where(m_cluster>=mass_selected[0])

	fo.close()


	plt.clf()

	plt.locator_params(nbins=16)

	if projection == "3D":
		ax1 = fig.add_subplot(111, projection='3d')
		ax1.set_zlabel('Z-Position (pc)')
		ax1.set_ylabel('Y-Position (pc)')
		ax1.set_zlim(-1*xy_box,xy_box)		
		ax1.set_ylim(-1*z_box,z_box)
		ax1.scatter(x,z, y, marker='.',color='black', s=2, alpha=1.)
		ax1.scatter(x_cluster[np.where(m_cluster>=mass_selected[0])],z_cluster[np.where(m_cluster>=mass_selected[0])], y_cluster[np.where(m_cluster>=mass_selected[0])], marker='*',color='red', s=20)	
					
	elif projection=="2D":
	
		ax1 = fig.add_subplot(111)
		ax1.set_ylabel('Y-Position (pc)')		
		ax1.set_ylim(-1*xy_box,xy_box)
		ax1.scatter(x, y, marker='.',color='black', s=20, alpha=0.3)
		
		ax1.scatter(x_cluster[np.where(m_cluster>=mass_selected[0])], y_cluster[np.where(m_cluster>=mass_selected[0])], marker='*',color='red', s=20)	
	

	ax1.set_xlabel('X-Position (pc)')
	ax1.set_xlim(-1*xy_box,xy_box)



	textstring='Time: %.2f Myr' % (round(float(time),2))
	plt.title(textstring)
	plt.grid(True)
        plt.text(0.75*xy_box,0.95*xy_box,"qvir = " + qvir)
        plt.text(0.75*xy_box,0.9*xy_box,"fdim = " + fdim)
        plt.text(-1.1*xy_box,-1.1*xy_box,"fbin = 0%",fontsize=11)
        plt.text(0.94*xy_box,-1.1*xy_box,kval,fontsize=10)


	if projection =="3D":
		plt.savefig(str(ifname)+'.xyz.png', dpi=my_dpi, bbox_inches='tight')
	elif projection == "2D":
		plt.savefig(str(ifname)+'.xy.png', dpi=my_dpi, bbox_inches='tight')
        plt.close()
	del(mass_srt_by_m)
	del(mass_selected)
        return snapname

arglist = sys.argv
sim = ''
if(len(arglist) < 2):

	sim = raw_input("Enter file path: ")
	#sim = 'analysis_test/outputs/runinv_k01'
else:

	sim = str(arglist[1])

kval = sim.split("_")[1] #get k01, k02, etc

picture_count = int(0)
nfile=0
fdim = raw_input("fdim: ")
qvir = raw_input("qvir: ")

ncores = int(raw_input("Enter number of cores (Enter -1 to turn off parallisation, 0 to autodetect system cores): "))

if (ncores>0):

		ncorearray=np.arange(0,ncores)
	
		ppservers=()
		job_server=pp.Server(ncpus=ncores,ppservers=ppservers)
	
		print '---------------------------------------------'
	
		print 'Number of cores initialised', ncorearray
		print 'Starting Parallel processing with', job_server.get_ncpus(), 'cores active'
	
		print '---------------------------------------------'
		
		time.sleep(2)

elif (ncores==0):
		print 'Autodetecting number of cores \n'

	
		ppservers=()
		job_server=pp.Server(ncpus='autodetect',ppservers=ppservers)
		ncorearray=np.arange(0,job_server.get_ncpus())	
		print '---------------------------------------------'
	
		print 'Number of cores initialised', ncorearray
		print 'Starting Parallel processing with', job_server.get_ncpus(), 'cores active'
	
		print '---------------------------------------------'
		
		time.sleep(2)
		
for fname in os.listdir(sim + '/snapshots/'):
    if 'snap' in fname and not 'tmp' in fname and not 'png' in fname:
        nfile += 1
print 'Number of files: ' + str(nfile)
picture_count = int(0)
ifname=[]
snapname=[]

for i in range(1,nfile+1):
	zeros=''
	for kx in range (0,(3 - int(np.floor(np.log10(i))))):
		zeros=zeros+str(0)
        ifname.append(str(sim + '/snapshots/' + '/snap'+ str(zeros) + str(i)))
	snapname.append('snap'+ str(zeros) + str(i))
	del zeros
	
	if ncores==-1:
   		generate_snapshot(ifname[i-1],snapname[i-1],projection,qvir,fdim,kval)

if ncores!=-1:   	
	irange = np.arange(0, i-1, 1)
	jobs = [(ix, job_server.submit(generate_snapshot,args=(ifname[ix],snapname[ix],projection,qvir,fdim,kval),depfuncs=(), modules=("matplotlib","numpy","mpl_toolkits"), globals=globals())) for ix in irange]
	for ix, job in jobs:
           px=job()
           del(px)
	   print snapname[ix]
	   #print "Snapshot", ifname[ix], "is", job()

os.system('avconv -y -r 24 -i ' + sim + '/snapshots/snap%04d.xy.png -s 1024x800 ' + sim + '/xy_' + kval + '.mp4')

for fn in os.listdir(sim + '/snapshots/'):
    if '.xy.png' in fn or '.xyz.png' in fn:
        os.system("rm -rf " + sim + '/snapshots/' + fn)
