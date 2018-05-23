library(ggplot2)
library(reshape2)
library(matrixStats) # for rowMedians

origin <- getwd()

fstr <- 'f16'
qstr <- 'q03'

fn3d <- 'allMSTs_lambar_3D'
proj <- 'xy' #currently doesn't do anything. Change to variable later.
fn2d <- 'allMSTs_lambar_xy'
cluster <- 'cluster_FoV5pc' # use 5 pc field of view

#'root' directory:
rootdir <- '/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0'

#select snapshots
#knum <- 'k01'
#snapshots <- c(1:220, 315:349)
knum <- 'k02'
snapshots <- c(490:520)
#knum <- 'k03'
#snapshots <- c(58, 77, 400, 425)
#knum <- 'k07'
#snapshots <- c(520:649)
#knum <- 'k10'
#snapshots <- c(25:42, 60:80, 120:135)

#input (data) and output (plots) directories:
inpath <- file.path(rootdir, paste0(fstr, qstr, '/analysis/runinv_', knum), cluster)
outpath <- file.path(rootdir, 'stats', paste0('pdf_', fstr, qstr, '_', knum))
#create output directory if it doesn't exist:
ifelse(!dir.exists(outpath), dir.create(outpath), FALSE)

nmsts <- 1000 #don't *need* this, it just helps with memory allocation. Might as well, since it's known

nsnaps <- length(readLines(file.path(inpath, 'CDFdata/allMSTs_lambar_3D.dat')))

all_lambdas3d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_3D.dat'),
                            nrows = nsnaps, row.names = 1, header = FALSE,
                            col.names = c("snap", sprintf("snap%04d", 1:nmsts)))

all_lambdas2d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_xy.dat'),
                            nrows = nsnaps, row.names = 1, header = FALSE,
                            col.names = c("snap", sprintf("snap%04d", 1:nmsts)))

#get data for each snapshot: nloop lambdas for each snapshot
# Order: [row, column] - data saved for rows specified in 'snapshots', all columns.
my_lambdas3d <- data.matrix(all_lambdas3d[snapshots, ])
my_lambdas2d <- data.matrix(all_lambdas2d[snapshots, ])

medians_3d <- data.frame(rowMedians(as.matrix(all_lambdas3d)))
medians_3d <- gsub("c|\\(|\\)|\n|,", "", medians_3d) #remove c, parentheses, \n, and commas from string
medians_2d <- data.frame(rowMedians(as.matrix(all_lambdas2d)))
medians_2d <- gsub("c|\\(|\\)|\n|,", "", medians_2d) #remove c, parentheses, \n, and commas from string


#kernel density estimates
#(see https://www.r-bloggers.com/the-density-function/ for good summary)
for (snap in 1:length(snapshots)){
  
  df_lambdas <- data.frame(my_lambdas3d[snap, ], my_lambdas2d[snap, ])
  colnames(df_lambdas)[1] <- '3D'; colnames(df_lambdas)[2] <- '2D'
  df_long <- melt(df_lambdas, variable.name = "Dimension", value.name = "Lambda") #convert data to long format
  
  pval <- wilcox.test(as.numeric(my_lambdas3d[snap, ]), as.numeric(my_lambdas2d[snap, ]), paired = FALSE)$p.value
  
  #set up output plot:
  outfn_pdf <- sprintf('snap%04dpdf_%02s.png', snapshots[snap], proj) # output file name
  png(filename = file.path(outpath, outfn_pdf), width = 500, height = 500)#, res=40)
  
  plot(
    ggplot(df_long, aes(x = Lambda, color = Dimension)) +
      geom_line(stat = "density") +
      scale_color_manual(values = c("black", "red")) +
      #lines at median 3D & 2D positions:
      geom_vline(xintercept = median(my_lambdas3d[snap, ]), colour = "black") +
      geom_vline(xintercept = median(my_lambdas2d[snap, ]), colour = "red", linetype = "dotted") +
      annotate("rect", xmin = median(my_lambdas3d[snap, ]) * 0.9, xmax = median(my_lambdas3d[snap, ]) * 1.1,
           ymin = -Inf, ymax = Inf, alpha = .2, fill = "blue") + #shaded region 10% around 3D med
      theme_minimal()  +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
      #annotate with p-value & snapshot number:
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = sprintf("p = %4.2e", pval)) +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 3.2, size = 3.6,
             label = ifelse(pval > 1.e-3, TRUE, FALSE)) +
      annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = sprintf("snap %04d", snapshots[snap]))
  )
  dev.off() #close plot
  
  #Plot lambda time-series showing this snapshot
  quartiles3d <- paste(quantile(df_lambdas[, 1]), collapse = ' ') #concaternate 'quantile' vector to one string
  quartiles2d <- paste(quantile(df_lambdas[, 2]), collapse = ' ')
  sixths3d <- paste(quantile(df_lambdas[, 1], probs = seq(0, 1, 1/6)), collapse = ' ') #concaternate 'quantile' vector to one string
  sixths2d <- paste(quantile(df_lambdas[, 2], probs = seq(0, 1, 1/6)), collapse = ' ')
  system(paste('python PDFs_lambda.py', outpath, knum, cluster, snapshots[snap], #model variables
               '--medians3d', medians_3d,
               '--medians2d', medians_2d,
               '--quartiles3d', quartiles3d,
               '--quartiles2d', quartiles2d,
               '--sixths3d', sixths3d,
               '--sixths2d', sixths2d
         ), wait=FALSE)
  
  #get y-values at given x-value at various points to find any snaps that don't overlap
  #(to check whether 2D & 3D 'agree' for paper purposes)
  # Get coordinates of density function for 3D & 2D lambdas
  dens3D <- density(df_lambdas[, 1])
  dens2D <- density(df_lambdas[, 2])
  dens3D_y <- approxfun(dens3D$x, dens3D$y)
  dens2D_y <- approxfun(dens2D$x, dens2D$y)
  
  #demonstrate method to find point on density function using median as example:
  hist(df_lambdas[, 1], freq = F)
  lines(dens3D, col="red", lwd=2)
  abline(v=median(df_lambdas[, 1]), lty=2)
  points(median(df_lambdas[, 1]), dens3D_y(median(df_lambdas[, 1])),
         cex=1.2, pch=20, col="blue")
  points(min(df_lambdas[, 1]), dens3D_y(min(df_lambdas[, 1])),
         cex=1.2, pch=20, col="blue")
  points(max(df_lambdas[, 1]), dens3D_y(max(df_lambdas[, 1])),
         cex=1.2, pch=20, col="blue")
  #Use 1/6 & 5/6 (0.17, 0.83) rather than quartiles:
  points(quantile(df_lambdas[, 1], 0.17), dens3D_y(quantile(df_lambdas[, 1], 0.17)),
         cex=1.2, pch=20, col="blue")
  points(quantile(df_lambdas[, 1], 0.83), dens3D_y(quantile(df_lambdas[, 1], 0.83)),
         cex=1.2, pch=20, col="blue")
  #draw 2D density function:
  lines(dens2D, col="black", lwd=2)
  #dev.off()
  
  #compare y-values of 2D & 3D density plots at max, min, median, & quartiles
  denspoints_3D <- c(dens3D_y(min(df_lambdas[, 1])),
                     dens3D_y(quantile(df_lambdas[, 1], .25)),
                     dens3D_y(median(df_lambdas[, 1])),
                     dens3D_y(quantile(df_lambdas[, 1], .75)),
                     dens3D_y(max(df_lambdas[, 1]))
  )#keep df_lambdas[, 1] for 2D check as x-values need to be the same
  denspoints_2D <- c(dens2D_y(min(df_lambdas[, 1])),
                     dens2D_y(quantile(df_lambdas[, 1], .25)),
                     dens2D_y(median(df_lambdas[, 1])),
                     dens2D_y(quantile(df_lambdas[, 1], .75)),
                     dens2D_y(max(df_lambdas[, 1]))
  )
}

#'ecdf': Empirical Cumulative Distribution Function
#for (snap in 1:length(snapshots)) {
#  CDF1 <- ecdf(my_lambdas3d[snap, ])
#  CDF2 <- ecdf(my_lambdas2d[snap, ])
#  outfn_cdf <- sprintf('lambdas%02s_snap%04d.png', proj, snapshots[snap]) # output file name
#  #set up output plot:
#  png(filename = file.path(outpath, outfn_cdf), width = 500, height = 500)#, res=40)
#  plot(CDF1, do.points = FALSE) #plot CDF1
#  lines(CDF2, col='red', do.points = FALSE) #overplot CDF2
#  legend("topleft", legend = c("3D", "xy"), col = c('black', 'red'), lty = c(1, 1))
#  dev.off() #close plot
#}