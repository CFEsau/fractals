library(scatterplot3d)
library(animation)
source('/local/cfe/backed_up_on_astro3/fractals/graphics/plots/3D/addgrids3d.r')

nkvals <- 10
for (k in 1:nkvals){
  knum <- sprintf('k%02d',k)
  
  fbin <- 'fbinary1p0'
  fdim <- 'f16'
  qvir <- 'q05'
  
  masterdir <- '/local/cfe/backed_up_on_astro3/fractals/r1p0'
  outdir <- file.path(masterdir,fbin,paste0(fdim,qvir),'analysis')
  snapdir <- file.path(outdir,paste0('runinv_',knum),'snapshots')
  setwd(snapdir)
  
  nmst <- 10
  axlim <- 6  #I'm not sure why the plot doesn't like '5'/certain numbers...
  fovlim <- 5 #Field of view limit in pc
  
  #Set delay between frames when replaying
  #ani.options(interval=0.03,loop=1)
  #saveVideo({
    
    #find number of snapshots
    nsnaps <- length(list.files(pattern="^snap"))
    #nsnaps <- 100
    #loop over all snapshots
    for (i in 1:nsnaps){
      
      #set input and output filenames
      infn <- sprintf('snap%04d',i)
      outfn <- sprintf('3D/3Dsnap%04d.png',i)
      #set up output plot
      #png(filename=outfn)
      
      #read data and save to vectors
      snapdata <- read.table(infn,col.names=c('IDdat','IDstar','tstar',
                                              'mstar','rx','ry','rz','vx','vy','vz'))
      rx <- snapdata$rx; ry <- snapdata$ry; rz <- snapdata$rz
      mstar <- snapdata$mstar
      
      #set up data frame with coordinates and masses ordered by decreasing mass
      star_df <- data.frame(rx,ry,rz,mstar)[order(-mstar),]
      
      #Define centre of cluster to be at mean x, y, z coordinates
      #so the cluster doesn't drift out of plotting space
      xmean <- mean(rx); ymean <- mean(ry); zmean <- mean(rz)
      
      #Get subset of data containing stars within FoV (centred around mean (x,y,z))
      plotstars_df <- subset(star_df, star_df[,1]<fovlim+xmean &
                               star_df[,1]>-fovlim+xmean &
                               star_df[,2]<fovlim+ymean & star_df[,2]>-fovlim+ymean &
                               star_df[,3]<fovlim+zmean & star_df[,3]>-fovlim+zmean)
      
      rx <- plotstars_df[,1]; ry <- plotstars_df[,2]; rz <- plotstars_df[,3]
      mstar <- plotstars_df[,4]
      
      #Make 3D scatter plot of star positions
      par(xpd=NA) #allows text outside of plotting region (see ?par)
      
      #1: Create empty grid
      spo <- scatterplot3d(star_df[,1:3], pch = "", grid=FALSE, box=FALSE,
                  xlim=c(-axlim,axlim),ylim=c(-axlim,axlim),zlim=c(-axlim,axlim),
                  xlab="x (pc)", ylab="",zlab="z (pc)")
      
      time <- (i/nsnaps)*10
      title(main = sprintf("Time = %.2f Myr",time),line=-2,adj=0)
      #add 'snap####' label, using 'par("usr")' vector to set position:
      #("usr": A vector of the form c(x1, x2, y1, y2) giving the
      # extremes of the user coordinates of the plotting region.)
      text(par("usr")[1]+0.8,par("usr")[4]*0.8,infn)
      
      #2: Add 3D grid lines
      addgrids3d(x=c(-axlim,axlim), y=c(-axlim,axlim), z=c(-axlim,axlim),
                 grid = c("xy","xz","yz"),
                 col.grid = "grey", lty.grid=par("lty"))
      
      #3: Add points
      spo$points3d(rx,ry,rz,type="p",pch=20)
      #generate points for projection on xy plane ('seq(...)' generates list of -axlim)
      spo$points3d(rx,ry,seq(-axlim,-axlim,length.out=length(plotstars_df[,1])),
                   type="p",pch=".")
      
      #add the nmst most massive stars as red dots
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
#  },video.name=paste0(outdir,"/movies/",knum,"_",fovlim,"pc.mp4"))
} #end of kvals