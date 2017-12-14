library(scatterplot3d)
library(animation)
source('/media/claire/Elements/Work/graphics/Rcode/addgrids3d.r')

knum='k03'
qvir='qvir0p3'
fdim='fdim1p6'
modelpath <- file.path(fdim,qvir,"outputs",paste0('runinv_',knum),'snapshots')
modeldir <- paste0('/media/claire/Elements/Work/fbin0p0/',modelpath)
setwd(modeldir)
#https://stackoverflow.com/questions/29680328/use-variable-in-file-path-in-r

nmst <- 10
axlim <- 6
fovlim=5 #Field of view limit in pc

#Set delay between frames when replaying
ani.options(interval=0.03,loop=1)
saveVideo({
  
  #find number of snapshots
  nsnaps <- length(list.files(pattern="^snap"))
  #nsnaps <- 100
  #loop over all snapshots
  for (i in 1:nsnaps){
    #set input and output filenames
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
    #set up output plot
    #png(filename=outfn)
    #read data and save to vectors
    snapdata <- read.table(infn)
    IDdat <- snapdata$V1
    IDstar<- snapdata$V2
    tstar <- snapdata$V3
    mstar <- snapdata$V4
    rx <-snapdata$V5
    ry <-snapdata$V6
    rz<-snapdata$V7
    vx <-snapdata$V8
    vy<-snapdata$V9
    vz<-snapdata$V10
    
    #Define centre of cluster to be at mean x, y, z coordinates
    xmean <- mean(rx)
    ymean <- mean(ry)
    zmean <- mean(rz)
    
    #set up data frame with coordinates and masses ordered by decreasing mass
    star_df <- data.frame(rx,ry,rz,mstar)[order(-mstar),]
    
    #Get subset of data containing stars within FoV (centred around mean (x,y,z))
    plotstars_df <- subset(star_df, star_df[,1]<fovlim+xmean & star_df[,1]>-fovlim+xmean & 
                           star_df[,2]<fovlim+ymean & star_df[,2]>-fovlim+ymean & 
                           star_df[,3]<fovlim+zmean & star_df[,3]>-fovlim+zmean)
    #(https://stackoverflow.com/questions/9648030/r-xlim-ylim-and-zlim-not-working-for-rgl-plot3d)
    rx <- plotstars_df[,1]
    ry <- plotstars_df[,2]
    rz <- plotstars_df[,3]
    mstar <- plotstars_df[,4]
    
    #Make 3D scatter plot of star positions
    #1: Create empty grid
    par(xpd=NA)
    spo <- scatterplot3d(star_df[,1:3], pch = "", grid=FALSE, box=FALSE,
                xlim=c(-axlim,axlim),ylim=c(-axlim,axlim),
                zlim=c(-axlim,axlim),xlab="x (pc)", ylab="",zlab="z (pc)")
    time <- (i/nsnaps)*10
    title(main = sprintf("Time = %.2f Myr",time),line=-2,adj=0)
    text(par("usr")[1]+0.8,par("usr")[4]*0.8,infn)
    #2: Add 3D grid lines
    addgrids3d(x=c(-axlim,axlim), y=c(-axlim,axlim),
               z=c(-axlim,axlim),grid = c("xy","xz","yz"),
               col.grid = "grey", lty.grid=par("lty"))
    #3: Add points
    spo$points3d(rx,ry,rz,type="p",pch=20)
    #generate points for projection on xy plane ('seq(...)' generates list of -axlim)
    spo$points3d(rx,ry,seq(-axlim,-axlim,length.out=length(plotstars_df[,1])),
                           type="p",pch=".")
    #add the nmst most massive stars as red dots
    #with(plotstars_df[plotstars_df$mstar >= m_sorted[nmst],],
    spo$points3d(rx[1:nmst],ry[1:nmst],rz[1:nmst],type="h",pch=19,cex=1,col='red')
    
    #Add text
    text(axlim*1.3, axlim*1.08,labels=knum,pos=2,cex=0.95)
    text(axlim*1.3, axlim*1.01,labels=qvir,pos=2,cex=0.95)
    text(axlim*1.3, axlim*0.93,labels=fdim,pos=2,cex=0.95)
    text(par("usr")[1],par("usr")[3]-0.6,labels="fbin=0%",cex=0.95)
    
    #for y axis label parallel with y axis:
    dims <- par("usr")
      ylab_x <- dims[1]+ 0.9*diff(dims[1:2])
      ylab_y <- dims[3]+ 0.1*diff(dims[3:4])
      text(ylab_x,ylab_y,"y (pc)",srt=45)
    #dev.off()
  }
},video.name="animatesnaps.mp4")