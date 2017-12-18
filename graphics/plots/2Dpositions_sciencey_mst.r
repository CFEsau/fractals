#!/usr/bin/Rscript
#www.statmethods.net/advgraphs/layout.html
#2D plots at given projection
library(animation)

nkvals <- 10
clustype <- 'cluster_FoV5pc' # one of _all, _FoV#pc, _r#rhalf

for (k in 1:nkvals){
  knum <- sprintf('k%02d',k)
  
  fbin='fbin0p0'
  fdim='f16'
  qvir='q03'
  
  masterdir <- '~/Documents/Work'
  outdir <- file.path(masterdir,fbin,paste0(fdim,qvir))
  kdir <- file.path(outdir,'outputs',paste0('runinv_',knum))
  snapdir <- file.path(kdir,'snapshots')
  clusterdir <- file.path(kdir,clustype)
  #paste0('/media/claire/Elements/Work/',fbin,fdim,qvir)
  #modeldir <- paste0('/media/claire/Elements/Work/fbin0p0/',fdim,qvir,'/outputs')
  #snapdir <- file.path(modeldir,paste0('runinv_',knum),'snapshots')
  setwd(snapdir)
  
  axlim <- 5.0
  fovlim <- 5.0 #Field of view limit in pc
  
  #find number of snapshots
  nsnaps <- length(list.files(pattern="^snap"))
  
  #loop over all projections
  for (p in 1:3){
    proj <- if (p==1) 'xy' else if (p==2) 'yz' else if (p==3) 'xz'
    
    ifelse(!dir.exists(file.path(snapdir, proj)),
           dir.create(file.path(snapdir, proj)), FALSE)
    
    #Movie: Set delay between frames when replaying
    #ani.options(interval=0.02,loop=1)
    #saveVideo({
    
    #loop over all snapshots
    for (i in 1:nsnaps){
      #set input and output filenames
      infn <- sprintf('snap%04d',i)
      
      #read data and save to vectors
      snapdata <- read.table(infn)
      mstar_all <- snapdata$V4
      rx_all <- snapdata$V5
      ry_all <- snapdata$V6
      rz_all <- snapdata$V7
      
      #Define centre of cluster to be at mean x, y, z coordinates
      xmean <- mean(rx_all)
      ymean <- mean(ry_all)
      zmean <- mean(rz_all)
      
      #set up data frame with coordinates and masses ordered by decreasing mass
      star_df <- data.frame(rx_all,ry_all,rz_all,mstar_all)[order(-mstar_all),]
      
      outfn <- sprintf('%02ssnap%04d.png',proj,i)
      #set up output plot
      #png(filename=file.path(proj,outfn),width = 500, height = 500)#,res=40)
      
      
      #Get subset of data containing stars within FoV (centred around mean (x,y,z))
      plotstars_df <- if (p==1) subset(star_df, star_df[,1]<fovlim+xmean & 
                                         star_df[,1]>-fovlim+xmean &
                                         star_df[,2]<fovlim+ymean &
                                         star_df[,2]>-fovlim+ymean) else if(
                                           p==2) subset(star_df, star_df[,2]<fovlim+ymean &
                                                          star_df[,2]>-fovlim+ymean &
                                                          star_df[,3]<fovlim+zmean &
                                                          star_df[,3]>-fovlim+zmean) else if(
                                                            p==3) subset(star_df, star_df[,1]<fovlim+xmean &
                                                                           star_df[,1]>-fovlim+xmean &
                                                                           star_df[,3]<fovlim+zmean &
                                                                           star_df[,3]>-fovlim+zmean)
      
      rx <- plotstars_df[,1] - xmean
      ry <- plotstars_df[,2] - ymean
      rz <- plotstars_df[,3] - zmean
      mstar <- plotstars_df[,4]
      
      #get object star positions
      obj_infn <- paste0('../',clustype,'/lambda/coords/',infn,'_objpositions_',proj,'.dat')
      obj_df <- read.table(obj_infn)
      objrx <- obj_df[,1] - xmean
      objry <- obj_df[,2] - ymean
      objrz <- obj_df[,3] - zmean
      
      
      #get MST edge coordinates
      mst_infn <- paste0('../',clustype,'/lambda/coords/',infn,'_objconnections_',proj,'.dat')
      mst_df <- read.table(mst_infn)
      x0 <- mst_df$V1
      y0 <- mst_df$V2
      z0 <- mst_df$V3
      x1 <- mst_df$V4
      y1 <- mst_df$V5
      z1 <- mst_df$V6
      
      #Make scatter plot of star positions
      
      #Projection:
      xcoord <- as.name(paste0("r",substr(proj,1,1)))
      ycoord <- as.name(paste0("r",substr(proj,2,2)))
      xdat <- eval(parse(text = as.name(xcoord)))
      ydat <- eval(parse(text = as.name(ycoord)))
      #object stars:
      xobjcoord <- as.name(paste0("objr",substr(proj,1,1)))
      yobjcoord <- as.name(paste0("objr",substr(proj,2,2)))
      xobjdat <- eval(parse(text = as.name(xobjcoord)))
      yobjdat <- eval(parse(text = as.name(yobjcoord)))
      
      
      #Create empty grid
      if (length(xcoord)<1) {xdat <- 0; ydat <- 0; xobjdat <- 0; yobjdat <- 0} #Set values to 0 if no points in data frame
      
      #Add stars:
      plot(xdat, ydat, col='gray60', cex=0.5,pch=20,
           #axes = FALSE,
	         #ann = FALSE,
          xlim=c(-axlim,axlim), ylim=c(-axlim,axlim), xlab=paste0(as.name(substr(proj,1,1)),' (pc)'),
              ylab=paste0(as.name(substr(proj,2,2)),' (pc)'), cex.lab=1)
      points(xobjdat,yobjdat,col='darkred',cex=1.2,pch=20)
      
        #if (p==1) {
        #  segments(x0[j], y0[j], x1[j], y1[j])#, col = par("fg"), lty = par("lty"), xpd = FALSE)
        #} else if (p==2) {
        #  segments(y0[j], z0[j], y1[j], z1[j])#, col = par("fg"), lty = par("lty"), xpd = FALSE)
        #} else if (p==3) {
        #  segments(x0[j], z0[j], x1[j], z1[j])#, col = par("fg"), lty = par("lty"), xpd = FALSE)
        #}
      if (p==1) {
        segments(x0, y0, x1, y1)#, col = par("fg"), lty = par("lty"), xpd = FALSE)
      } else if (p==2) {
        segments(y0, z0, y1, z1)#, col = par("fg"), lty = par("lty"), xpd = FALSE)
      } else if (p==3) {
        segments(x0, z0, x1, z1)#, col = par("fg"), lty = par("lty"), xpd = FALSE)
      }
      
      mtext(sprintf('t = %.2f Myr',(i/nsnaps)*10),side=3,adj=1)#,col.main = "gray80",cex.main=0.9)
      
      #dev.off() #close plot
    }#end of snapshot loop
    #},video.name=paste0(outdir,"/movies/",knum,proj,"_",fovlim,"pc.mp4")) #end of video
  }#end of projection loop
}#end of knum loop
