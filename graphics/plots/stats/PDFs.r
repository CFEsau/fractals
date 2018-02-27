library(ggplot2)
library(reshape2)

origin <- getwd()

fn3d <- 'allMSTs_lambar_3D'
proj <- 'xy' #currently doesn't do anything. Change to variable later.
fn2d <- 'allMSTs_lambar_xy'
cluster <- 'cluster_FoV5pc'
outpath <- paste0('/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0/f16q03/analysis')

#select snapshots
#knum <- 'k01'
#snapshots <- c(76, 105, 135, 459)
#knum <- 'k02'
#snapshots <- c(200, 275, 400, 420)
#knum <- 'k03'
#snapshots <- c(58, 77, 400, 425)
knum <- 'k10'
snapshots <- c(25:42,60:80,120:135)

simpath <- file.path(outpath,paste0('runinv_',knum),cluster)
plotsdir <- paste0(outpath,'/plots/cdf_',knum,'_test') #outputs directory for plots
#create output directory if it doesn't exist:
ifelse(!dir.exists(plotsdir), dir.create(plotsdir), FALSE)

nmsts <- 1000
macro3d <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_3D.dat'), row.names=1, nrows=nmsts)
macro2d <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_xy.dat'), row.names=1, nrows=nmsts)

#get data for each snapshot: nloop lambdas for each snapshot
# Order: [row, column] - data saved for rows specified in 'snapshots', all columns.
lambdas3d <- data.matrix(macro3d[snapshots,])
lambdas2d <- data.matrix(macro2d[snapshots,])

#'ecdf': Empirical Cumulative Distribution Function
#for (snap in 1:length(snapshots)){
#  CDF1 <- ecdf(lambdas3d[snap,])
#  CDF2 <- ecdf(lambdas2d[snap,])
#  outfn <- sprintf('lambdas%02s_snap%04d.png',proj,snapshots[snap]) # output file name
#  #set up output plot:
#  png(filename=file.path(plotsdir,outfn),width = 500, height = 500)#,res=40)
#  plot(CDF1, do.points=FALSE) #plot CDF1
#  lines(CDF2, col='red', do.points=FALSE) #overplot CDF2
#  legend("topleft",legend=c("3D","xy"),col=c('black','red'),lty=c(1,1))
#  dev.off() #close plot
#}

#kernel density estimates
#(see https://www.r-bloggers.com/the-density-function/ for good summary)
for (snap in 1:length(snapshots)){
  
  df <- data.frame(lambdas3d[snap,],lambdas2d[snap,])
  colnames(df)[1] <- '3D'; colnames(df)[2] <- '2D'
  df_long <- melt(df) #convert data to long format
  
  pval <- wilcox.test(as.numeric(lambdas3d[snap,]),as.numeric(lambdas2d[snap,]),paired=FALSE)$p.value

  plot(
    ggplot(df_long, aes(x=value,color=variable)) +
      geom_line(stat="density") +
      scale_color_manual(values=c("black","red")) +
      geom_vline(xintercept=median(lambdas3d[snap,]),colour="black") + #line at median 3D point
      geom_vline(xintercept=median(lambdas2d[snap,]),colour="red",linetype="dotted") + #line at median 3D point
      annotate("rect",xmin=median(lambdas3d[snap,])*0.9,xmax=median(lambdas3d[snap,])*1.1,
           ymin=-Inf,ymax=Inf,alpha=.2,fill="blue") + #shaded region 10% around 3D med
      theme_minimal()  +
      scale_x_continuous(breaks=scales::pretty_breaks(n=6)) +
      #annotate with p-value & snapshot number:
      annotate("text",x=Inf,y=Inf,hjust=1,vjust=1,label=sprintf("p = %4.2e",pval)) +
      annotate("text",x=Inf,y=Inf,hjust=1,vjust=3.2,size=3.6,
             label=ifelse(pval>1.e-3,TRUE,FALSE)) +
      annotate("text",x=-Inf,y=Inf,hjust=0,vjust=1,label=sprintf("snap %04d",snapshots[snap]))
  )
  
  ##density data (nloop lambdas) taken from each row of data table (1st PDF from 1st snapshot, etc).
  #dens3d <- density(lambdas3d[snap,])
  #dens2d <- density(lambdas2d[snap,])
  #outfn2 <- sprintf('pdf%02s_snap%04d.png',proj,snapshots[snap]) # output file name
  ##set up output plot:
  #png(filename=file.path(plotsdir,outfn2),width = 500, height = 500)#,res=40)
  #ylim <- max(range(dens3d$y, dens2d$y))
  #plot(dens3d,ylim=c(0, ylim),main=snapshots[snap])
  #lines(dens2d, col='red')
  #legend("topright",legend=c("3D","xy"),col=c('black','red'),lty=c(1,1))
  
  ##calculate the p-value for the distribution & annotate
  #pval <- wilcox.test(as.numeric(lambdas3d[snap,]), as.numeric(lambdas2d[snap,]),paired=FALSE)$p.value
  #usr <- par("usr")
  #text(usr[1]*1.02,usr[4]*0.99,labels=sprintf("p = %4.2e",pval),adj=c(0,1))
  #text(usr[1]*1.02,usr[4]*0.92,labels=ifelse(pval>1.e-3,TRUE,FALSE),adj=c(0,0))
  
  ##add vertical line at median lambda
  #med3d <- median(lambdas3d[snap,])
  #abline(v=med3d,col="blue")
  ##add shaded region to +- 10% median
  #rect(xleft=med3d*0.9,ybottom=usr[3],xright=med3d*1.1,ytop=usr[4],col="lightblue",lty=0)
  
  #dev.off() #close plot
}



#ks.test(lambdas3d[3,],lambdasxy[3,])

#wilcox.test(lambdas3d[snap,],lambdasxy[snap,])
#pfromU <- wilcox.test(lambdas3d[snap,],lambdasxy[snap,])

#t.test(lambdas3d[snap,],lambdasxy[snap,])
#pfromt <- t.test(lambdas3d[snap,],lambdasxy[snap,])