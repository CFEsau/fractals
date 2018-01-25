library(ggplot2)

origin <- getwd()

fn3d <- 'allMSTs_lambar_3D'
proj <- 'xy' #currently doesn't do anything. Change to variable later.
fn2d <- 'allMSTs_lambar_xy'
cluster <- 'cluster_FoV5pc'
outpath <- paste0('/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary1p0/f16q05/analysis')

#select snapshots
#knum <- 'k01'
#snapshots <- c(76, 105, 135, 459)
knum <- 'k02'
snapshots <- c(200, 275, 400, 420)
#knum <- 'k03'
#snapshots <- c(58, 77, 400, 425)

simpath <- file.path(outpath,paste0('runinv_',knum),cluster)
plotsdir <- paste0(outpath,'/plots/cdf_',knum) #outputs directory for plots
#create output directory if it doesn't exist:
ifelse(!dir.exists(plotsdir), dir.create(plotsdir), FALSE)

nmsts <- 1000
macro3d <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_3D.dat'), row.names=1, nrows=nmsts)
macro2d <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_xy.dat'), row.names=1, nrows=nmsts)

#get data for each snapshot
# Order: [row, column] - data saved for rows specified in 'snapshots', all columns.
lambdas3d <- data.matrix(macro3d[snapshots,])
lambdas2d <- data.matrix(macro2d[snapshots,])

#'ecdf': Empirical Cumulative Distribution Function
for (snap in 1:length(snapshots)){
  CDF1 <- ecdf(lambdas3d[snap,])
  CDF2 <- ecdf(lambdas2d[snap,])
  outfn <- sprintf('lambdas%02s_snap%04d.png',proj,snapshots[snap]) # output file name
  #set up output plot:
  png(filename=file.path(plotsdir,outfn),width = 500, height = 500)#,res=40)
  plot(CDF1, do.points=FALSE) #plot CDF1
  lines(CDF2, col='red', do.points=FALSE) #overplot CDF2
  legend("topleft",legend=c("3D","xy"),col=c('black','red'),lty=c(1,1))
  dev.off() #close plot
}

#kernel density estimates
#(see https://www.r-bloggers.com/the-density-function/ for good summary)
for (snap in 1:length(snapshots)){
  dens3d <- density(lambdas3d[snap,])
  densxy <- density(lambdasxy[snap,])
  outfn2 <- sprintf('pdf%02s_snap%04d.png',proj,snapshots[snap]) # output file name
  #set up output plot:
  png(filename=file.path(plotsdir,outfn2),width = 500, height = 500)#,res=40)
  ylim <- max(range(dens3d$y, densxy$y))
  plot(dens3d,ylim=c(0, ylim))
  lines(densxy, col='red')
  legend("topright",legend=c("3D","xy"),col=c('black','red'),lty=c(1,1))
  dev.off() #close plot
}



#ks.test(lambdas3d[3,],lambdasxy[3,])

wilcox.test(lambdas3d[snap,],lambdasxy[snap,])
pfromU <- wilcox.test(lambdas3d[snap,],lambdasxy[snap,])

t.test(lambdas3d[snap,],lambdasxy[snap,])
pfromt <- t.test(lambdas3d[snap,],lambdasxy[snap,])