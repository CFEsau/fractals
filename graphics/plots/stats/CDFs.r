library(ggplot2)

origin <- getwd()

fdim <- 'fdim1p6'
qvir <- 'qvir0p3'

fn3d <- 'allMSTs_lambar_3D'
fnxy <- 'allMSTs_lambar_xy'
cluster <- 'cluster_FoV5pc'
outpath <- paste0('../../../fbin0p0/f16q03/outputs')


#knum <- 'k01'
#snapshots <- c(76, 105, 135, 459)
knum <- 'k02'
snapshots <- c(200, 275, 400, 420)
#knum <- 'k03'
#snapshots <- c(58, 77, 400, 425)

listnum <- 4

simpath <- file.path(outpath,paste0('runinv_',knum),cluster)
plotsdir <- paste0(outpath,'/plots/cdf_',knum) #outputs directory for plots

nmsts <- 1000
macro3d <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_3D.dat'), row.names=1, nrows=nmsts)
macroxy <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_xy.dat'), row.names=1, nrows=nmsts)

#get data for each snapshot ('snapshots' object). Row #, all columns.
# Order: [row, column]
lambdas3d <- data.matrix(macro3d[snapshots,])
lambdasxy <- data.matrix(macroxy[snapshots,])
#plot(ecdf(lambdas3d[listnum,]))
CDF1 <- ecdf(lambdas3d[listnum,])
CDF2 <- ecdf(lambdasxy[listnum,])
plot(CDF1)
lines(CDF2, col='red')

#CDF1 <- ecdf(macro3d[snapshots,]) # data for selected snapshots
#CDF2 <- ecdf(macroxy[snapshots,])
#ks.test(lambdas3d[3,],lambdasxy[3,])
wilcox.test(lambdas3d[listnum,],lambdasxy[listnum,])
pfromU <- wilcox.test(lambdas3d[listnum,],lambdasxy[listnum,])

t.test(lambdas3d[listnum,],lambdasxy[listnum,])
pfromt <- t.test(lambdas3d[listnum,],lambdasxy[listnum,])

dens3d <- density(lambdas3d[listnum,])
densxy <- density(lambdasxy[listnum,])
ylim <- max(range(dens3d$y, densxy$y))
plot(dens3d,ylim=c(0, ylim))
lines(densxy, col='red')
