#!/usr/bin/Rscript
#www.statmethods.net/advgraphs/layout.html
#2D plots at given projection
library(ggplot2)
library(animation)
library(Hmisc) #for minor tick marks
library("plotrix") #for draw.circle


#function for drawing a circle on plot:
#(if using plot instead of ggplot could just use draw.circle, but 'plot' doesn't
#like my axis limits & end up with non-uniform gaps between circle x,y & axes)
#could also use geom_circle in ggforce but can't install ggforce package - 
# need dev version of udunits2 package (need to do sudo apt-get install libudunits2-dev)
circleFun <- function(center = c(0,0), r = 0.5, npoints = 100){
  tt <- seq(0,2*pi,length.out=npoints)
  xx <- center[1] + r*cos(tt)
  yy <- center[2] + r*sin(tt)
  return(data.frame(xcir=xx, ycir=yy))
}



nkvals <- 10
clustype <- "cluster_FoV5pc" # one of _all, _FoV#pc, _r#rhalf

fbin <- "fbinary0p0"
fvals <- c(1.6, 2.0, 2.6, 3.0); fstr <- c("f16", "f20", "f26", "f30")
qvals <- c(0.3, 0.5); qstr <- c("q03", "q05")

masterdir <- "/local/cfe/backed_up_on_astro3/fractals/r1p0"

#use larger axis limits if using all stars in cluster
axlim <- c(5.0,5.0) #c(10,10)  #1st & current limit (not x & y!) Used for dynamicallim. 
                    #limits not working... 4.7 actually gives ~5
dynamicallim <- FALSE  #increase/decrease axis limits with time? Initial limits given by axlim
fovlim <- 5.0         #Field of view limit in pc
usefovlim <- ifelse(clustype!="cluster_all",TRUE,FALSE) #don't restrict massive stars selection
                                                        #to field of view if using 'cluster_all'

theme_set(theme_bw() + theme(panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())) #pre-set the bw theme for plots


for (f in 1:length(fvals)) {
  fdim <- fstr[f]
  
  for (q in 1:length(qvals)) {
    qvir <- qstr[q]
    outdir <- file.path(masterdir,fbin,paste0(fdim,qvir),'analysis')
    message(file.path(fbin,paste0(fdim,qvir),'analysis')) #print model directory (no prefix)
    
    for (k in 1:nkvals){
      knum <- sprintf('k%02d',k)
      message(sprintf("k = %d",k))
      
      kdir <- file.path(outdir,paste0('runinv_',knum))
      snapdir <- file.path(kdir,"snapshots")
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
        ani.options(interval=0.02,loop=1)
        saveVideo({
          
          #loop over all snapshots
          for (i in 1:nsnaps){
            infn <- sprintf('snap%04d',i) #input file name
            
            #read data and save to vectors
            snapdata <- read.table(infn)
            
            #set up data frame with coordinates and masses ordered by decreasing mass
            star_df <- data.frame("rx_all"=snapdata$V5,"ry_all"=snapdata$V6,"rz_all"=snapdata$V7,
                                   "mstar_all"=snapdata$V4)[order(-snapdata$V4),]
            
            #Define centre of cluster to be at mean x, y, z coordinates
            xmean <- mean(star_df$rx_all)
            ymean <- mean(star_df$ry_all)
            zmean <- mean(star_df$rz_all)
            
            outfn <- sprintf('%02ssnap%04d.png',proj,i) #output file name
            #set up output plot
            #png(filename=file.path(proj,outfn),width = 500, height = 500)#,res=40)
            
            #Get subset of data containing stars within FoV (centred around mean (x,y,z))
            ifelse(usefovlim,
                   plotstars_df <- if (
                     p==1) subset(star_df, sqrt( (star_df[,1]-xmean)^2 + (star_df[,2]-ymean)^2 ) < fovlim) else if(
                       p==2) subset(star_df, sqrt( (star_df[,2]-ymean)^2 + (star_df[,3]-zmean)^2 ) < fovlim) else if(
                         p==3) subset(star_df, sqrt( (star_df[,1]-xmean)^2 + (star_df[,3]-zmean)^2 ) < fovlim),
                   plotstars_df <- star_df
            )
            
            colnames(plotstars_df) <- c("rx","ry","rz","mstar")
            plotstars_df$rx <- plotstars_df$rx - xmean
            plotstars_df$ry <- plotstars_df$ry - ymean
            plotstars_df$rz <- plotstars_df$rz - zmean
            
            #dynamical axis limits: maximum absolute position coordinate
            if (dynamicallim){
              maxcoord <- if (proj=='xy') max(abs(c(rx,ry))) else if (
                proj=='yz') max(abs(c(ry,rz))) else if (
                  proj=='xz') max(abs(c(rx,rz)))
              axlim[2] <- ifelse(maxcoord>axlim[1],maxcoord,axlim[1]) #axlim 1 is initial axis limit (i.e. minimum limit)
            }
            
            #get object star positions
            obj_infn <- paste0('../',clustype,'/lambda/coords/',infn,'_objpositions_',proj,'.dat')
            obj_df <- read.table(obj_infn)
            colnames(obj_df) <- c("objrx","objry","objrz")
            obj_df$objrx <- obj_df$objrx - xmean
            obj_df$objry <- obj_df$objry - ymean
            obj_df$objrz <- obj_df$objrz - zmean
            
            
            #Make scatter plot of star positions
            #Plot stars:
            gg <- ggplot(plotstars_df,aes(
              x=eval(parse(text=
                             as.name(paste0("plotstars_df$r",substr(proj,1,1))))), #get x coord from 1st field of 'proj' (rx/ry/rz)
              y=eval(parse(text=
                             as.name(paste0("plotstars_df$r",substr(proj,2,2))))) #get y coord from second field of 'proj'
            )) +
              geom_point(colour = "gray60",size=0.5) +
              
              # overplot 'object' stars:
              geom_point(data = obj_df, aes(
                x=eval(parse(text = 
                               as.name(paste0("obj_df$objr",substr(proj,1,1))))),
                y=eval(parse(text = 
                          as.name(paste0("obj_df$objr",substr(proj,2,2)))))
                ), 
                colour = "darkred") +
              
              # define tick mark positions:
              scale_x_continuous(breaks=scales::pretty_breaks(n=6)) +
              scale_y_continuous(breaks=scales::pretty_breaks(n=6)) +
              labs(subtitle=as.name(sprintf('t = %.2f Myr',(i/nsnaps)*10)), # t= ... Myr
                   x=paste0(as.name(substr(proj,1,1)),' (pc)'),
                   y=paste0(as.name(substr(proj,2,2)),' (pc)'))
            
            # add circle at FoV limit
            if (usefovlim){
              gg <- gg + geom_path(data=circleFun(r=5),aes(x=xcir,y=ycir),
                                   colour = "gray80",linetype="dashed")
              }
            
            plot(gg)
            
            #dev.off() #close plot
            
          }#end of snapshot loop
        },video.name=paste0(outdir,"/movies/",knum,proj,"_",fovlim,"pc_test2.mp4")) #end of video
      }#end of projection loop
    }#end of knum loop
  }#end of qvals loop
}#end of fvals loop
