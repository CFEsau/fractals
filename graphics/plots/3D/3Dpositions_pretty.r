library(scatterplot3d)
#library(animation)
source('/local/cfe/backed_up_on_astro3/fractals/graphics/plots/3D/addgrids3d.r')

nkvals <- 10
for (k in 1:nkvals){
  if (k<10) {
    knum=paste0('k0',k)}
  else {
    knum=paste0('k',k)}
  
  fbin='fbinary1p0'
  fdim='f16'
  qvir='q05'

  masterdir <- '/local/cfe/backed_up_on_astro3/fractals/r1p0'
  outdir <- file.path(masterdir,fbin,paste0(fdim,qvir),'analysis')
  kdir <- file.path(outdir,paste0('runinv_',knum))
  snapdir <- file.path(kdir,'snapshots')
  setwd(snapdir)

  nmst <- 10
  axlim <- 5.0
  fovlim <- 5.0 #Field of view limit in pc

#Set delay between frames when replaying
#ani.options(interval=0.02,loop=1)
#saveVideo({
  
  #find number of snapshots
  nsnaps <- length(list.files(pattern="^snap"))
  
  #loop over all snapshots
  for (i in 1:nsnaps){
    
    #set input and output filenames
    infn <- sprintf('snap%04d',i)
    outfn <- sprintf('3D/3Dsnap%04d.png',i)
    #set up output plot
    #png(filename=outfn)
    
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
    
    
    #Get subset of data containing stars within FoV (centred around mean (x,y,z))
    plotstars_df <- subset(star_df, star_df[,1]<fovlim+xmean & star_df[,1]>-fovlim+xmean & 
                           star_df[,2]<fovlim+ymean & star_df[,2]>-fovlim+ymean & 
                           star_df[,3]<fovlim+zmean & star_df[,3]>-fovlim+zmean)
    #(https://stackoverflow.com/questions/9648030/r-xlim-ylim-and-zlim-not-working-for-rgl-plot3d)
    
    rx <- plotstars_df[,1] - xmean
    ry <- plotstars_df[,2] - ymean
    rz <- plotstars_df[,3] - zmean
    mstar <- plotstars_df[,4]
    colstar <- apply(plotstars_df,1,colourFunction)
    #pchstar <- apply(plotstars_df,1,pchFunction)
    cexstar <- apply(plotstars_df,1,cexFunction)
    
 
    #Make 3D scatter plot of star positions
    par(bg = 'black') # black background
    
    #colour palette:
    #spectype<-colorRampPalette(c("red","orange","yellow","#FFFFCC","blue")) # red, orange, beige, blue
    
    
    #1: Create empty grid
    par(xpd=NA)
    if (length(rx)>0)  #Only add points if there is data within the FoV. Otherwise leave as empty grid.
      spo <- scatterplot3d(rx,y=ry,z=rz,pch = "", grid=FALSE, box=FALSE,
                  xlim=c(-axlim,axlim), ylim=c(-axlim,axlim), zlim=c(-axlim,axlim))
    
    #2: Add 3D grid lines
    addgrids3d(x=c(-axlim,axlim), y=c(-axlim,axlim),
               z=c(-axlim,axlim),grid = c("xy","xz","yz"),
               col.grid = "gray50", lty.grid=par("lty"))
    
    
    #3: Add points
    if(length(rx)>0)
      #generate points for projection on xy plane ('seq(...)' generates list of values at 'axlim')
      spo$points3d(rx,ry,seq(-axlim,-axlim,length.out=length(rz)),         
                   type="p",pch=".", col='gray30')
      # yz plane:
      spo$points3d(seq(-axlim,-axlim,length.out=length(rx)),ry,rz,
                   type="p",pch=".", col='gray30')
      # xz plane:
      spo$points3d(rx,seq(axlim,axlim,length.out=length(ry)),rz,
                   type="p",pch=".", col='gray30')
      
      
      #Add most massive stars on projections as +s with different colour
      #spo$points3d(rx[1:nmst],ry[1:nmst],seq(-axlim,-axlim,length.out=nmst),type="p",pch='+',col='purple')
      #spo$points3d(rx[1:nmst],seq(axlim,axlim,length.out=nmst),rz[1:nmst],type="p",pch='+',col='purple')
      #spo$points3d(seq(-axlim,-axlim,length.out=nmst),ry[1:nmst],rz[1:nmst],type="p",pch='+',col='purple')
      
      #Add stars
      #spo$points3d(rx,ry,rz,type="p",pch=20, col=spectype(mstar))
      spo$points3d(rx,ry,rz, type='p', pch=20, cex=cexstar, col=colstar)
      #col=ifelse(((abs(X)>1.65 & abs(Y)>1.65)),"red", "black")
      
    time <- (i/nsnaps)*10
    title(main = sprintf("Time = %.2f Myr",time),line=-2,adj=0.1,col.main = "gray80", cex.main = 0.9)
      
    #dev.off()
  }
#},video.name=paste0(outdir,"/movies/",knum,"_",fovlim,"pc.mp4"))
}
