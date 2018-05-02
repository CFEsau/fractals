circleFun <- function(center = c(0,0), r = 0.5, npoints = 100){
  tt <- seq(0,2*pi,length.out=npoints)
  xx <- center[1] + r*cos(tt)
  yy <- center[2] + r*sin(tt)
  return(data.frame(xcir=xx, ycir=yy))
}

plotstars <- function(x, y, starcol='gray60', xobj=NULL, yobj=NULL, objcol='darkred', limit, xvar, yvar, snap, nsnaps,
                      plotmst=FALSE,xmst0=NULL,xmst1=NULL,ymst0=NULL,ymst1=NULL,
                      fovlim=5,plotfov=FALSE){
  
  plot(xdat, ydat, col=starcol, cex=0.5, pch=20,
       xlim=c(-limit,limit), ylim=c(-limit,limit), xlab=paste0(xvar,' (pc)'), ylab=paste0(yvar,' (pc)'), cex.lab=1)
  points(xobj,yobj,col=objcol,cex=1.2,pch=20)
  minor.tick(nx=2,ny=2,tick.ratio=0.4)
  mtext(sprintf('t = %.2f Myr',(i/nsnaps)*10),side=3,adj=1)#,col.main = "gray80",cex.main=0.9)
  
  if (plotfov){
    #add circle for FoV limit (could use circleFun function as commented out. Keeping original for posterity. Speed difference?)
    #circle <- circleFun(r=fovlim)
    #points(circle$xcir,circle$ycir)
    curve(sqrt(25-x^2),-fovlim,fovlim,n=200,add=TRUE,type="l",lty=2,col='gray80')
    curve(-sqrt(25-x^2),-fovlim,fovlim,n=200,add=TRUE,type="l",lty=2,col='gray80')
  }
  
  #use 'segments' to draw a line between pairs of points
  if (plotmst) {
    segments(xmst0, ymst0, xmst1, ymst1)#, col = par("fg"), lty = par("lty"), xpd = FALSE)
  }
  
}



ggplotstars <- function(stardf = plotstars_df, objdf=obj_df, xvar, yvar, plotfov=TRUE,fovlim=5,
                        plotmst=FALSE,edgecoords=NULL){
  
  gg <- ggplot(stardf,aes_string(x=as.name(paste0("r",xvar)),
                                 y=as.name(paste0("r",yvar))))+ #rx/ry/rz from data frame
  
    geom_point(colour = "gray60",size=0.5) +
    
    # overplot 'object' stars:
    geom_point(data = obj_df, aes_string(x=as.name(paste0("r",xvar)),y=as.name(paste0("r",yvar))),colour = "darkred") +
    
    # define tick mark positions:
    scale_x_continuous(breaks=scales::pretty_breaks(n=6)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=6)) +
    labs(subtitle=as.name(sprintf('t = %.2f Myr',(i/nsnaps)*10)), # t= ... Myr
         x=paste0(as.name(substr(proj,1,1)),' (pc)'),
         y=paste0(as.name(substr(proj,2,2)),' (pc)'))
  
  # add mst edges
  if (plotmst){
    gg <- gg + geom_segment(data=plotedges,aes(x=x0,y=y0,xend=x1,yend=y1))
  }
  
  # add circle at FoV limit
  if (plotfov){
    gg <- gg + geom_path(data=circleFun(r=fovlim),aes(x=xcir,y=ycir),
                         colour = "gray80",linetype="dashed")
  }
  
  
  
  plot(gg)
}