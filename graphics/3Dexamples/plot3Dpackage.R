library("plot3D") #Load plot3D package
#http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization

setwd('/media/claire/Elements/Work/fbin0p0/f16q03/outputs/runinv_k01/snapshots') #Data directory
nmst <- 10 #Number of stars in mst
axlim <- 0.5 #Axis limits (pc)

nsnaps <- length(list.files(pattern="^snap")) #Count number of snapshots in data directory

#for (i in 1:nsnaps){ #Loop over all snapshots
i=1
  #set input and output filenames:
  if (i<10) {
    infn=paste('snap000',i,sep='')
    outfn=paste('snap000',i,'.png',sep='')}
  if (i<100 && i>=10) {
    infn=paste('snap00',i,sep='')
    outfn=paste('snap00',i,'.png',sep='')}
  if (i<1000 && i>=100) {
    infn=paste('snap0',i,sep='')
    outfn=paste('snap0',i,'.png',sep='')}
  if (i<10000 && i>=1000) {
    infn=paste('snap',i,sep='')
    outfn=paste('snap',i,'.png',sep='')}
  
  snapdata <- read.table(infn) #Read in data (V# are default column names)
  mstar <- snapdata$V4
  rx <-snapdata$V5
  ry <-snapdata$V6
  rz<-snapdata$V7
  
  #star_df=data.frame(rx,ry,rz,mstar) #Set up data fram with coordinates and masses
  star_df <- data.frame(rx,ry,rz,mstar)[order(-mstar),]
  scatter3D(rx, ry, rz, bty = "b2",colkey = FALSE, clab = c("rz (pc)"),
            xlim=c(-axlim,axlim),ylim=c(-axlim,axlim),zlim=c(-axlim,axlim),
            ticktype = "detailed",phi=10,theta=30)
#}