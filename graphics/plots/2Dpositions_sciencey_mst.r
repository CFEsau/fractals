#!/usr/bin/Rscript
#www.statmethods.net/advgraphs/layout.html
#2D plots at given projection
library(animation)
library(Hmisc) #for minor tick marks

nkvals <- 10
clustype <- 'cluster_all' # one of _all, _FoV#pc, _r#rhalf


fbin='fbinary0p0'
fvals <- c(1.6, 2.0, 2.6, 3.0); fstr <- c("f16", "f20", "f26", "f30")
qvals <- c(0.3, 0.5); qstr <- c("q03", "q05")

masterdir <- '/local/cfe/backed_up_on_astro3/fractals/r1p0'

#use larger axis limits if using all stars in cluster
axlim <- c(10,10)  #1st & current limit, not x & y! Used for dynamicallim. #limits not working... 4.7 actually gives ~5
dynamicallim <- TRUE  #increase/decrease axis limits with time. Initial limits given by axlim
fovlim <- 5.0         #Field of view limit in pc
usefovlim <- ifelse(clustype!="cluster_all",TRUE,FALSE) #don't restrict massive stars selection
                                                        #to field of view if using 'cluster_all'

for (f in 1:length(fvals)) {

  for (q in 1:length(qvals)) {
    
    fdim <- fstr[f]; qvir <- qstr[q]
    message(file.path(fbin,paste0(fdim,qvir),'analysis'))
    
    for (k in 1:nkvals){
      
      knum <- sprintf('k%02d',k)
      message(sprintf("k = %d",k))
      
      outdir <- file.path(masterdir,fbin,paste0(fdim,qvir),'analysis')
      kdir <- file.path(outdir,paste0('runinv_',knum))
      snapdir <- file.path(kdir,'snapshots')
      clusterdir <- file.path(kdir,clustype)
      setwd(snapdir)
      
      #find number of snapshots
      nsnaps <- length(list.files(pattern="^snap"))
      
      #loop over all projections
      for (p in 1:3){
        
        proj <- if (p==1) 'xy' else if (p==2) 'yz' else if (p==3) 'xz'
        message(sprintf("\t%s...",proj))
        
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
            
            outfn <- sprintf('mst%02ssnap%04d.png',proj,i)
            #set up output plot
            png(filename=file.path(proj,outfn),width = 500, height = 500)#,res=40)
            
            ifelse(usefovlim,
                   plotstars_df <- if (
                     p==1) subset(star_df, sqrt( (star_df[,1]-xmean)^2 + (star_df[,2]-ymean)^2 ) < fovlim) else if(
                       p==2) subset(star_df, sqrt( (star_df[,2]-ymean)^2 + (star_df[,3]-zmean)^2 ) < fovlim) else if(
                         p==3) subset(star_df, sqrt( (star_df[,1]-xmean)^2 + (star_df[,3]-zmean)^2 ) < fovlim),
                   plotstars_df <- star_df
            )
            
            #centre stars around mean:
            rx <- plotstars_df[,1] - xmean
            ry <- plotstars_df[,2] - ymean
            rz <- plotstars_df[,3] - zmean
            mstar <- plotstars_df[,4]
            
            #axis limits: maximum absolute position coordinate
            if (dynamicallim){
              maxcoord <- if (proj=='xy') max(abs(c(rx,ry))) else if (
                proj=='yz') max(abs(c(ry,rz))) else if (
                  proj=='xz') max(abs(c(rx,rz)))
              axlim[2] <- ifelse(maxcoord>axlim[1],maxcoord,axlim[1]) #axlim 1 is initial axis limit (i.e. minimum limit)
            }
            
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
                 xlim=c(-axlim[2],axlim[2]), ylim=c(-axlim[2],axlim[2]), xlab=paste0(as.name(substr(proj,1,1)),' (pc)'),
                 ylab=paste0(as.name(substr(proj,2,2)),' (pc)'), cex.lab=1)
            points(xobjdat,yobjdat,col='darkred',cex=1.2,pch=20)
            minor.tick(nx=2,ny=2,tick.ratio=0.4)
            if (usefovlim){
              #add circle for FoV limit (have to do as 2 curves, because R...)
              curve(sqrt(25-x^2),-5,5,n=200,add=TRUE,type="l",lty=2,col='gray80')
              curve(-sqrt(25-x^2),-5,5,n=200,add=TRUE,type="l",lty=2,col='gray80')
            }
            
            #use 'segments' to draw a line between pairs of points
            if (p==1) {
              segments(x0, y0, x1, y1)#, col = par("fg"), lty = par("lty"), xpd = FALSE)
            } else if (p==2) {
              segments(y0, z0, y1, z1)#, col = par("fg"), lty = par("lty"), xpd = FALSE)
            } else if (p==3) {
              segments(x0, z0, x1, z1)#, col = par("fg"), lty = par("lty"), xpd = FALSE)
            }
            
            mtext(sprintf('t = %.2f Myr',(i/nsnaps)*10),side=3,adj=1)#,col.main = "gray80",cex.main=0.9)
            
            dev.off() #close plot
            
          }#end of snapshot loop
        #},video.name=paste0(outdir,"/movies/",knum,proj,"_",fovlim,"pc_new.mp4")) #end of video
      }#end of projection loop
    }#end of knum loop
  }#end of qvals loop
}#end of fvals loop
