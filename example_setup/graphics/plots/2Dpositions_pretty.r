#!/usr/bin/Rscript
#www.statmethods.net/advgraphs/layout.html
#2D plots at given projection
library(animation)

nkvals <- 10

for (k in 1:nkvals){
  knum <- sprintf('k%02d',k)
  
  fbin='fbin0p0'
  fdim='f16'
  qvir='q05'
  
  masterdir <- '/media/claire/Elements/Work'
  outdir <- file.path(masterdir,fbin,paste0(fdim,qvir),'outputs')
  kdir <- file.path(outdir,paste0('runinv_',knum))
  snapdir <- file.path(kdir,'snapshots')
  setwd(snapdir)
  
  nmst <- 10
  axlim <- 5.0
  fovlim <- 5.0 #Field of view limit in pc
  
  #find number of snapshots
  nsnaps <- length(list.files(pattern="^snap"))
  
  #loop over all projections
  for (p in 1:3){
    proj <- if (p==1) 'xy' else if (p==2) 'yz' else if (p==3) 'xz'
    
    ifelse(!dir.exists(file.path(snapdir, proj)),
           dir.create(file.path(snapdir, proj)), FALSE)
    
    #Set delay between frames when replaying
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
    
    # add another column with star colours:
    colourFunction <- function(x){
      m <- x[4]
      colour <- if(m<0.15) 'red3' else if (
        m>=0.15 && m<0.3) 'orangered' else if (
          m>=0.3 && m<0.9) 'darkorange3' else if (
            m>=0.9 && m<1.3) 'gold2' else if (
              m>=1.3 && m<1.6) 'lightyellow' else if (
                m>=1.6 && m<2.4) 'lightcyan' else if (
                  m>=2.4 && m<15) 'skyblue1' else if (
                    m>=15) 'deepskyblue'
      return(colour)
    }
    
    
    #add another column with cex scaling factor
    cexFunction <- function(x){
      
      m <- x[4]
      cexfac <- if(m<1.66) 0.5*(1.06*m^0.945) else if ( #0.856 at m=1.66
        m>=1.66) 0.5*(1.37*m^0.45) #WRONG! but need smaller high-m symbols
          #m>=1.66) 0.5*(1.33*m^0.555)
      return(cexfac)
    }
    
    #add another column with cex scaling factor
    #cexFunction <- function(x){
    #  m <- x[4]
    #  cexfac <- if(m<0.9) 0.5 else if (
    #    m>=0.9 && m<1.6) 0.7 else if (
    #      m>=1.6 && m<2.4) 0.9 else if (
    #        m>=2.4 && m<4) 1 else if (
    #          m>=4 && m<8) 1.2 else if (
    #            m>=8 && m<15) 1.5 else if (
    #              m>=15) 1.8
    #  return(cexfac)
    #}
    
    ##set up layout for multiple plots
    #layout(matrix(c(1,1,2,3),2,2,byrow=TRUE),widths=c(1,1),heights=c(1,1))
    ##layout.show(n=3)
    #par(mgp=c(0.2,0,0),mar=c(1.2,1.2,1.2,1.2),pty="s")
    
      
      outfn <- sprintf('pretty%02ssnap%04d_fdim.png',proj,i) # output file name
      #set up output plot
      #png(filename=outfn,width = 720, height = 720)
      
      #Get subset of data containing stars within FoV (centred around mean (x,y,z))
      #plotstars_df <- subset(star_df, star_df[,1]<fovlim+xmean & star_df[,1]>-fovlim+xmean & 
      #                       star_df[,2]<fovlim+ymean & star_df[,2]>-fovlim+ymean & 
      #                       star_df[,3]<fovlim+zmean & star_df[,3]>-fovlim+zmean)
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
      colstar <- apply(plotstars_df,1,colourFunction)
      #pchstar <- apply(plotstars_df,1,pchFunction)
      cexstar <- apply(plotstars_df,1,cexFunction)
      
      #Make scatter plot of star positions
      
      #Projection:
      xcoord <- as.name(paste0("r",substr(proj,1,1)))
      ycoord <- as.name(paste0("r",substr(proj,2,2)))
      xdat <- eval(parse(text = as.name(xcoord)))
      ydat <- eval(parse(text = as.name(ycoord)))
      
      #Create empty grid
      par(bg = 'black')
      if (length(xcoord)<1) {xdat <- 0; ydat <- 0} #Set xdat, ydat to 0 if there are no points in the data frame
      
      #Add stars:
      plot(xdat, ydat, col=colstar, cex=cexstar, pch=20,
           #axes = FALSE,
	         #ann = FALSE,
          xlim=c(-axlim,axlim), ylim=c(-axlim,axlim), col.lab="white",
               xlab=as.name(substr(proj,1,1)), ylab=as.name(substr(proj,2,2)),
    	           cex.lab=1.4)
          axis(1, col = "white", lwd.ticks=0)
      axis(2, col = "white", lwd.ticks=0)
    
      time <- (i/nsnaps)*10
      title(main = sprintf("Time = %.2f Myr",time),line=-26,adj=0.9,col.main = "gray80", cex.main = 0.8)
      #dev.off()
      
  }#end of snapshot loop
    #},video.name=paste0(outdir,"/movies/",knum,proj,"_",fovlim,"pc.mp4")) #end of video
  }#end of projection loop
}#end of knum loop

