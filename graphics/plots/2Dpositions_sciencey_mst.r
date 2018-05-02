#!/usr/bin/Rscript
#www.statmethods.net/advgraphs/layout.html
#2D plots at given projection

#inputs from 2Dpositions.sh:
args <- commandArgs()
#print(args)
f <- as.numeric(args[6]); q <- as.numeric(args[7])#; k <- as.numeric(args[8])

library(animation)
library(Hmisc) #for minor tick marks
source('/local/cfe/backed_up_on_astro3/Github/
       fractals/graphics/plots/plots2Dfn.r')


clustype <- "FoV5pc" # one of all, FoV#pc, r#rhalf
#(check obj star selection if using rhalf... may only work for 'all' & 'FoV')
plotmst = TRUE

axlimStart <- 5       #axis limits for initial plot
dynamicallim <- FALSE  #increase/decrease axis limits with time

fovlim <- 5        #Field of view limit in pc
usefovlim <- ifelse(clustype!="all",TRUE,FALSE)
#(don't restrict massive star selection to fov lims if using 'cluster_all')

ifelse(dynamicallim, plottype <- 'plot', plottype <- 'ggplot')
#Haven't tested ggplot with dynamical axis limits so don't allow for now

fbin <- "fbinary0p0"
fvals <- c(1.6, 2.0, 2.6, 3.0); fstr <- c("f16", "f20", "f26", "f30")
qvals <- c(0.3, 0.5); qstr <- c("q03", "q05")
proj <- c('xy','yz','xz')

nkvals <- 10

masterpath <- "/local/cfe/backed_up_on_astro3/fractals/r1p0"

# Set up global theme for ggplots:
theme_set(theme_bw() + #dark-on-light ggplot2 theme
            theme(panel.grid.major = element_blank(), #don't draw major
                  panel.grid.minor = element_blank(), # or minor grid lines
                  text = element_text(size=14) #pre-set plot text size
            )) #end of theme setup

#for (f in 1:length(fvals)) {
  fdim <- fstr[f]
  message("fdim: ", fdim)
  
  #for (q in 1:length(qvals)) {
    qvir <- qstr[q]
    
    modeldir <- file.path(masterpath,fbin,paste0(fdim,qvir),'analysis')
    message(file.path(fbin,paste0(fdim,qvir),'analysis')) #model directory
    
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
        psystime1 <- Sys.time() #note system time
        #message(psystime1)
        message(sprintf("\t%s...",proj[p]))
        
        #Get x & y coordinates for plot from 'proj' str
        xcoord <- substr(proj[p], 1, 1) #syntax: substr(string, from, to)
        ycoord <- substr(proj[p], 2, 2)
        
        ifelse(dynamicallim,
               outdir <- file.path(snapdir, sprintf(
                 '%s_dynamical_%s', proj[p], clustype)), #dynamical limits
               outdir <- file.path(snapdir, sprintf(
                 '%s_%s', proj[p], clustype)))  #fixed limits
        message(sprintf("\toutput directory: %s", basename(outdir)))
        
        ifelse(!dir.exists(outdir), dir.create(outdir), FALSE)
        
        #Movie: Set delay between frames when replaying
        #ani.options(interval=0.02,loop=1)
        #saveVideo({
          
          #loop over all snapshots
          for (i in 1:nsnaps){
            axlim <- axlimStart #new snapshot so reset axis limits
            
            infn <- sprintf('snap%04d', i) #input file name
            
            #read data and save to vectors
            snapdata <- read.table(infn)
            
            #data frame with coordinates & masses ordered by decreasing mass
            starDF <- data.frame("rx" = snapdata$V5,
                                 "ry" = snapdata$V6,
                                 "rz" = snapdata$V7,
                                 "m" = snapdata$V4)[order(-snapdata$V4), ]
            rm(snapdata) #remove snapdata to free up memory
            
            #Overwrite starDF with only the projections in 'proj':
            starDF <- starDF[c(paste0("r", xcoord), paste0("r", ycoord), "m")]
            
            #Define centre of cluster to be at mean position coordinates
            means <- c(mean(starDF[, 1]), mean(starDF[, 2]))
            
            #Subset of data containing stars within FoV centred around mean
            inFoV <- sqrt(
              (starDF[,1] - means[1])^2  + (starDF[,2] - means[2])^2 ) < fovlim
            ifelse(usefovlim,
                   plotstarsDF <- starDF[inFoV,], #plot star when FoV is TRUE
                   plotstarsDF <- starDF #plot all stars
            )
            rm(inFoV, starDF)
            
            #centre stars around mean:
            plotstarsDF[,1] <- plotstarsDF[,1] - means[1]
            plotstarsDF[,2] <- plotstarsDF[,2] - means[2]
            
            #dynamical axis limits: maximum absolute position coordinate
            if (dynamicallim){
              maxcoord <- max(abs(c(plotstarsDF[,1],plotstarsDF[,2])))
              axlim <- ifelse(maxcoord>axlim,maxcoord,axlim)
            }
            
            #get object star positions
            objFile <- paste0('../cluster_',clustype,'/lambda/coords/',
                              infn,'_objpositions_',proj[p],'.dat')
            objDF <- read.table(objFile)
            colnames(objDF) <- c("rx","ry","rz")
            
            #again only keep the dimensions in 'proj':
            objDF <- objDF[c(paste0("r",xcoord),paste0("r",ycoord))]
            
            #centre stars around mean:
            objDF[,1] <- objDF[,1] - means[1]
            objDF[,2] <- objDF[,2] - means[2]
            
            
            #Make scatter plot of star positions
            
            outfn <- sprintf('mst_%02ssnap%04d.png',proj[p],i)
            
            #set up output plot
            png(filename=file.path(outdir,outfn),
                width = 500, height = 500)#,res=40)
            
            if (plotmst) {
              #get MST edge coordinates
              mst_infn <- paste0('../cluster_',clustype,'/lambda/coords/',infn,
                                 '_objconnections_',proj[p],'.dat')
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
            #Projection: x coord from 1st field of 'proj' & y coord from 2nd
            xdat <- eval(parse(text = as.name(paste0("plotstarsDF$r",xcoord))))
            ydat <- eval(parse(text = as.name(paste0("plotstarsDF$r",ycoord))))
            #object stars:
            xobjdat <- eval(parse(text = as.name(paste0("objDF$r",xcoord))))
            yobjdat <- eval(parse(text = as.name(paste0("objDF$r",ycoord))))
            
            
            
              #Plot stars using 'plot':
              plot.stars(x=xdat, y=ydat, xobj=xobjdat, yobj=yobjdat,limit=axlim,
                        xvar=xcoord, yvar=ycoord, snap=i, nsnaps=nsnaps,
                        #plotmst=TRUE, xmst0=plotedges$x0,xmst1=plotedges$x1,
                        #ymst0=plotedges$y0,ymst1=plotedges$y1,
                        plotfov=usefovlim)
              
            } else if (plottype=='ggplot'){
              
              ggplot.stars(stardf=plotstarsDF, objdf=objDF,
                           xvar=xcoord, yvar=ycoord, plotfov=TRUE,fovlim=5,
                           plotmst=TRUE, edgecoords=plotedges)
              
            }
            
            dev.off() #close plot
            rm(plotstarsDF)
            
          }#end of snapshot loop
        #},video.name=paste0(modeldir,"/movies/",
                            #knum,proj[p],"_",fovlim,"pc_mst.mp4"))
        psystime2 <- Sys.time()
        #message(psystime2)
        message("Time elapsed in projection loop: ", psystime2-psystime1)
      }#end of projection loop
      ksystime2 <- Sys.time() #end timer
      #message(ksystime2)
      message("\n\tTime elapsed in k loop: ", ksystime2-ksystime1)
    }#end of knum loop
#  }#end of qvals loop
#}#end of fvals loop
