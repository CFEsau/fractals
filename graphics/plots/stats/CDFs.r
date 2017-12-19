library(ggplot2)

origin <- getwd()

fn3d <- 'allMSTs_lambar_3D'
fnxy <- 'allMSTs_lambar_xy'
cluster <- 'cluster_FoV5pc'
outpath <- paste0('../../../r0p5/fbin0p0/f16q03/analysis')

#select snapshots
#knum <- 'k01'
#snapshots <- c(76, 105, 135, 459)
knum <- 'k02'
snapshots <- c(200, 275, 400, 420)
#knum <- 'k03'
#snapshots <- c(58, 77, 400, 425)

simpath <- file.path(outpath,paste0('runinv_',knum),cluster)
plotsdir <- paste0(outpath,'/plots/cdf_',knum) #outputs directory for plots

nmsts <- 1000
macro3d <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_3D.dat'), row.names=1, nrows=nmsts)
macroxy <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_xy.dat'), row.names=1, nrows=nmsts)

#get data for each snapshot
# Order: [row, column] - data saved for rows specified in 'snapshots', all columns.
lambdas3d <- data.matrix(macro3d[snapshots,])
lambdasxy <- data.matrix(macroxy[snapshots,])

#'ecdf': Empirical Cumulative Distribution Function
for (snap in 1:length(snapshots)){
  CDF1 <- ecdf(lambdas3d[snap,])
  CDF2 <- ecdf(lambdasxy[snap,])
  plot(CDF1, do.points=FALSE) #plot CDF1
  lines(CDF2, col='red', do.points=FALSE) #overplot CDF2
}

#ks.test(lambdas3d[3,],lambdasxy[3,])

wilcox.test(lambdas3d[snap,],lambdasxy[snap,])
pfromU <- wilcox.test(lambdas3d[snap,],lambdasxy[snap,])

t.test(lambdas3d[snap,],lambdasxy[snap,])
pfromt <- t.test(lambdas3d[snap,],lambdasxy[snap,])

dens3d <- density(lambdas3d[snap,])
densxy <- density(lambdasxy[snap,])
ylim <- max(range(dens3d$y, densxy$y))
plot(dens3d,ylim=c(0, ylim))
lines(densxy, col='red')
