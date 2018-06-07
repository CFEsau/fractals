library(scatterplot3d)
library(animation)
source('/local/cfe/backed_up_on_astro3/fractals/graphics/plots/3D/addgrids3d.r')


fbin <- "fbinary0p0"
fvals <- c(1.6, 2.0, 2.6, 3.0); fstr <- c("f16", "f20", "f26", "f30")
qvals <- c(0.3, 0.5); qstr <- c("q03", "q05")
nkvals <- 10
clustype <- "cluster_all" # one of _all, _FoV#pc, _r#rhalf

nmst <- 10 #number of stars in MST
axlim <- c(6,6)  #Minimum limit & current limit. #I'm not sure why the plot doesn't like '5'/certain numbers...
fovlim <- 5 #Field of view limit in pc
usefovlim <- ifelse(clustype!="cluster_all",TRUE,FALSE) #don't restrict massive stars selection
                                                        #to field of view if using 'cluster_all'
dynamicallim <- TRUE #TRUE if axis limits change between snapshots, FALSE if static (stays as axlim)



masterdir <- "/local/cfe/backed_up_on_astro3/fractals/r1p0"

for (f in 1:length(fvals)) {
  fdim <- fstr[f]
  
for (q in 1:length(qvals)) {
    qvir <- qstr[q]
    outdir <- file.path(masterdir,fbin,paste0(fdim,qvir),'analysis')
    message(file.path(fbin,paste0(fdim,qvir),'analysis')) #print model directory (no prefix)
    
for (k in 1:nkvals){
      knum <- sprintf('k%02d',k)
      message(sprintf("k = %d",k))
      
      snapdir <- file.path(outdir,paste0('runinv_',knum),'snapshots')
      setwd(snapdir)
      ifelse(!dir.exists('3D'), dir.create('3D'), FALSE)
      
      #Set delay between frames when replaying
      #ani.options(interval=0.03,loop=1)
      #saveVideo({
      
      #find number of snapshots
      nsnaps <- length(list.files(pattern="^snap"))
      
      #loop over all snapshots
      for (i in 1:nsnaps){
        #for (i in 860:995){
        
        #set input and output filenames
        infn <- sprintf('snap%04d',i)
        outfn <- sprintf('3D/3Dsnap%04d.png',i)
        #set up output plot
        png(filename=outfn)
        
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
        ifelse(usefovlim, #if using fovlim then use stars within FoV;
               plotstars_df <- subset(star_df, star_df[,1]<fovlim+xmean &
                                        star_df[,1]>-fovlim+xmean &
                                        star_df[,2]<fovlim+ymean & star_df[,2]>-fovlim+ymean &
                                        star_df[,3]<fovlim+zmean & star_df[,3]>-fovlim+zmean),
               plotstars_df <- star_df    #else use all stars
        )
        
        rx <- plotstars_df[,1]; ry <- plotstars_df[,2]; rz <- plotstars_df[,3]
        mstar <- plotstars_df[,4]
        
        #dynamical axis limits: maximum absolute position coordinate
        if (dynamicallim){
          maxcoord <- max(abs(c(rx,ry,rz)))
          #axlim[2] <- ifelse(maxcoord>axlim[1],maxcoord,axlim[1]) #axlim 1 is initial axis limit (i.e. minimum limit)
          #
          # Axlim doesn't work properly with 3D axes... instead of setting max(abs(rx,ry,rz)),
          # find which values axlim is happy with and set these equal manually
          #
          if(maxcoord < 6.0){
            axlim[2] <- 6.0
          } else if (maxcoord >= 6.0 & maxcoord < 10.0){
            axlim[2] <- 10.0
          } else if (maxcoord >= 10.0 & maxcoord < 15.0){
            axlim[2] <- 15.0
          } else if (maxcoord >= 15.0 & maxcoord < 20.0){
            axlim[2] <- 20.0
          } else if (maxcoord >= 20.0 & maxcoord < 30.0){
            axlim[2] <- 30.0
          } else if (maxcoord >= 30.0 & maxcoord < 40.0){
            axlim[2] <- 40.0
          } else if (maxcoord >= 40.0 & maxcoord < 60.0){
            axlim[2] <- 60.0
          } else if (maxcoord >= 60.0 & maxcoord < 100.0){
            axlim[2] <- 100.0
          } else if (maxcoord >= 100.0 & maxcoord < 150.0){
            axlim[2] <- 150.0
          } else if (maxcoord >= 150.0 & maxcoord < 200.0){
            axlim[2] <- 200.0
          }  else if (maxcoord >= 200.0) {
            message("WARNING: maxcoord > 200. Increase axis limits accordingly.")
          }
          
        }
        #message(i,' ',maxcoord,' ',axlim[2])
        #Make 3D scatter plot of star positions
        par(xpd=NA) #allows text outside of plotting region (see ?par)
        
        #1: Create empty grid
        spo <- scatterplot3d(star_df[,1:3], pch = "", grid=FALSE, box=FALSE,
                             xlim=c(-axlim[2],axlim[2]),ylim=c(-axlim[2],axlim[2]),zlim=c(-axlim[2],axlim[2]),
                             xlab="x (pc)", ylab="",zlab="z (pc)")
        
        time <- (i/nsnaps)*10
        title(main = sprintf("Time = %.2f Myr",time),line=-2,adj=0)
        text(par("usr")[1]+0.8,par("usr")[4]*0.8,infn) #add 'snap####' label
        #'par("usr")': vector to set position ("usr": A vector of the form c(x1, x2, y1, y2)
        # giving the extremes of the user coordinates of the plotting region.)
        
        #2: Add 3D grid lines
        addgrids3d(x=c(-axlim[2],axlim[2]), y=c(-axlim[2],axlim[2]), z=c(-axlim[2],axlim[2]),
                   grid = c("xy","xz","yz"),
                   col.grid = "grey", lty.grid=par("lty"))
        
        #3: Add points
        spo$points3d(rx,ry,rz,type="p",pch=20)
        #generate points for projection on xy plane ('seq(...)' generates list of -axlim)
        spo$points3d(rx,ry,seq(-axlim[2],-axlim[2],length.out=length(plotstars_df[,1])),
                     type="p",pch=".")
        
        #add the nmst most massive stars as red dots
        spo$points3d(rx[1:nmst],ry[1:nmst],rz[1:nmst],type="h",pch=19,cex=1,col='red')
        
        #Add text
        text(par("usr")[2]*1.05, par("usr")[4],labels=knum,pos=2,cex=0.95)
        text(par("usr")[2]*1.05, par("usr")[4]*0.94,labels=qvir,pos=2,cex=0.95)
        text(par("usr")[2]*1.05, par("usr")[4]*0.88,labels=fdim,pos=2,cex=0.95)
        text(par("usr")[1],par("usr")[3]-0.6,labels="fbin=0%",cex=0.95)
        
        #for y axis label parallel with y axis:
        dims <- par("usr")
        ylab_x <- dims[1]+ 0.9*diff(dims[1:2])
        ylab_y <- dims[3]+ 0.1*diff(dims[3:4])
        text(ylab_x,ylab_y,"y (pc)",srt=45)
        
        dev.off()
      }
      #  },video.name=paste0(outdir,"/movies/",knum,"_",fovlim,"pc.mp4"))
    } #end of knum loop
  }#end of qvals loop
}#end of fvals loop
