#!/usr/bin/Rscript
#www.statmethods.net/advgraphs/layout.html
#2D plots at given projection
library(animation)
library(Hmisc) #for minor tick marks
source('/local/cfe/backed_up_on_astro3/Github/fractals/graphics/plots/plots2Dfn.r')


clustype <- "FoV5pc" # one of all, FoV#pc, r#rhalf (check obj star selection if using rhalf... may only work for 'all' & 'FoV')
plotmst = TRUE

axlim_start <- 5       #axis limits for initial plot
dynamicallim <- FALSE  #increase/decrease axis limits with time

fovlim <- 5        #Field of view limit in pc
usefovlim <- ifelse(clustype!="all",TRUE,FALSE) #don't restrict massive stars selection to field of view if using 'cluster_all'

ifelse(dynamicallim, plottype <- 'plot', plottype <- 'ggplot') #either plot or ggplot.
                                  #Haven't tested ggplot with dynamical axis limits so don't allow for now

fbin <- "fbinary0p0"
fvals <- c(1.6, 2.0, 2.6, 3.0); fstr <- c("f16", "f20", "f26", "f30")
qvals <- c(0.3, 0.5); qstr <- c("q03", "q05")
proj <- c('xy','yz','xz')

nkvals <- 10

masterpath <- "/local/cfe/backed_up_on_astro3/fractals/r1p0"

# Set up global theme for ggplots:
theme_set(theme_bw() + #dark-on-light ggplot2 theme
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), #don't draw major & minor grid lines
                  text = element_text(size=14) #pre-set plot text size
            )) #end of theme setup

for (f in 1:length(fvals)) {
  fdim <- fstr[f]
  message("fdim: ", fdim)
  
  for (q in 1:length(qvals)) {
    qvir <- qstr[q]
    
    modeldir <- file.path(masterpath,fbin,paste0(fdim,qvir),'analysis')
    message(file.path(fbin,paste0(fdim,qvir),'analysis')) #print model directory (no prefix)
    
    for (k in 1:nkvals){
      
      #set a timer for each k:
      ksystime1 <- Sys.time()
      #message(ksystime1)
      
      knum <- sprintf('k%02d',k)
      message(sprintf("k = %d",k))
      
      snapdir <- file.path(modeldir,paste0('runinv_',knum),"snapshots")
      setwd(snapdir)
      
      #find number of snapshots
      nsnaps <- length(list.files(pattern="^snap"))
      
      #loop over all projections
      for (p in 1:3){
        psystime1 <- Sys.time() #note system time (used when finding why code slows down over iterations)
        #message(psystime1)
        message(sprintf("\t%s...",proj[p]))
        
        #Get x & y coordinates for plot from 'proj' str
        xcoord <- substr(proj[p],1,1) #syntax: substr(string, from, to)
        ycoord <- substr(proj[p],2,2)
        
        ifelse(dynamicallim,
               outdir <- file.path(snapdir, sprintf('%s_dynamical_%s',proj[p],clustype)), #dynamical limits
               outdir <- file.path(snapdir, sprintf('%s_%s',proj[p],clustype)))  #defined fixed limits
        message(sprintf("\toutput directory: %s",basename(outdir)))
        
        ifelse(!dir.exists(outdir), dir.create(outdir), FALSE)
        
        #Movie: Set delay between frames when replaying
        #ani.options(interval=0.02,loop=1)
        #saveVideo({
          
          #loop over all snapshots
          for (i in 1:nsnaps){
            axlim <- axlim_start #new snapshot, so reset axis limits (put here in case limits decrease)
            
            infn <- sprintf('snap%04d',i) #input file name
            
            #read data and save to vectors
            snapdata <- read.table(infn)
            
            #set up data frame with coordinates and masses ordered by decreasing mass
            star_df <- data.frame("rx"=snapdata$V5,"ry"=snapdata$V6,"rz"=snapdata$V7,
                                  "m"=snapdata$V4)[order(-snapdata$V4),]
            rm(snapdata) #remove snapdata to free up memory
            
            #Define centre of cluster to be at mean rx/ry/rz
            xmean <- mean(star_df$rx)
            ymean <- mean(star_df$ry)
            zmean <- mean(star_df$rz)
            
            #Get subset of data containing stars within FoV (centred around mean projection coordinates)
            ifelse(usefovlim,
                   plotstars_df <- if (
                     p==1) subset(star_df, sqrt( (star_df[,1]-xmean)^2 + (star_df[,2]-ymean)^2 ) < fovlim) else if(
                       p==2) subset(star_df, sqrt( (star_df[,2]-ymean)^2 + (star_df[,3]-zmean)^2 ) < fovlim) else if(
                         p==3) subset(star_df, sqrt( (star_df[,1]-xmean)^2 + (star_df[,3]-zmean)^2 ) < fovlim),
                   plotstars_df <- star_df
            )
            
            rm(star_df)
            
            colnames(plotstars_df) <- c("rx","ry","rz","m")
            #centre stars around mean:
            plotstars_df$rx <- plotstars_df$rx - xmean
            plotstars_df$ry <- plotstars_df$ry - ymean
            plotstars_df$rz <- plotstars_df$rz - zmean
            
            #dynamical axis limits: maximum absolute position coordinate
            if (dynamicallim){
              maxcoord <- if (proj[p]=='xy') max(abs(c(plotstars_df$rx,plotstars_df$ry))) else if (
                proj[p]=='yz') max(abs(c(plotstars_df$ry,plotstars_df$rz))) else if (
                  proj[p]=='xz') max(abs(c(plotstars_df$rx,plotstars_df$rz)))
              axlim <- ifelse(maxcoord>axlim,maxcoord,axlim)
              #axlim[2] <- ifelse(maxcoord>axlim[1],maxcoord,axlim[1]) #axlim 1 is initial axis limit (i.e. minimum limit)
            }
            
            #get object star positions
            obj_infn <- paste0('../cluster_',clustype,'/lambda/coords/',infn,'_objpositions_',proj[p],'.dat')
            obj_df <- read.table(obj_infn)
            colnames(obj_df) <- c("rx","ry","rz")
            obj_df$rx <- obj_df$rx - xmean
            obj_df$ry <- obj_df$ry - ymean
            obj_df$rz <- obj_df$rz - zmean
            
            
            #Make scatter plot of star positions
            
            outfn <- sprintf('mst_%02ssnap%04d.png',proj[p],i) #output file name
            #set up output plot
            png(filename=file.path(outdir,outfn),width = 500, height = 500)#,res=40)
            
            if (plotmst) {
              #get MST edge coordinates
              mst_infn <- paste0('../cluster_',clustype,'/lambda/coords/',infn,'_objconnections_',proj[p],'.dat')
              mst_df <- read.table(mst_infn)
              colnames(mst_df) <- c("x0","y0","z0","x1","y1","z1")
              
              plotedges <- data.frame(
                "x0"=eval(parse(text = as.name(paste0("mst_df$",xcoord,"0")))),
                "y0"=eval(parse(text = as.name(paste0("mst_df$",ycoord,"0")))),
                "x1"=eval(parse(text = as.name(paste0("mst_df$",xcoord,"1")))),
                "y1"=eval(parse(text = as.name(paste0("mst_df$",ycoord,"1"))))
              )
            }
            
            if (plottype=='plot'){
            #Projection:
            xdat <- eval(parse(text = as.name(paste0("plotstars_df$r",xcoord)))) #get x coord from 1st field of 'proj' (xy/xz/yz)
            ydat <- eval(parse(text = as.name(paste0("plotstars_df$r",ycoord)))) #get y coord from second field of 'proj'
            #object stars:
            xobjdat <- eval(parse(text = as.name(paste0("obj_df$r",xcoord))))
            yobjdat <- eval(parse(text = as.name(paste0("obj_df$r",ycoord))))
            
            
            
              #Plot stars using 'plot':
              plotstars(x=xdat, y=ydat, xobj=xobjdat, yobj=yobjdat,limit=axlim, xvar=xcoord, yvar=ycoord, snap=i, nsnaps=nsnaps,
                        #plotmst=TRUE, xmst0=plotedges$x0,xmst1=plotedges$x1,ymst0=plotedges$y0,ymst1=plotedges$y1,
                        plotfov=usefovlim)
              
            } else if (plottype=='ggplot'){
              
              ggplotstars(stardf=plotstars_df, objdf=obj_df, xvar=xcoord, yvar=ycoord, plotfov=TRUE,fovlim=5,
                          plotmst=TRUE, edgecoords=plotedges)
              
            }
            
            dev.off() #close plot
            rm(plotstars_df)
            
          }#end of snapshot loop
        #},video.name=paste0(modeldir,"/movies/",knum,proj[p],"_",fovlim,"pc_mst.mp4")) #end of video
        psystime2 <- Sys.time()
        #message(psystime2)
        message("Time elapsed in projection loop: ", psystime2-psystime1)
      }#end of projection loop
      ksystime2 <- Sys.time() #end timer
      #message(ksystime2)
      message("\n\tTime elapsed in k loop: ", ksystime2-ksystime1)
    }#end of knum loop
  }#end of qvals loop
}#end of fvals loop
