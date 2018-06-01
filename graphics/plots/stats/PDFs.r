#Produce probability density functions for user-defined snapshots

library(ggplot2)
library(reshape2)
library(matrixStats) # for rowMedians
library(dplyr) # for 'lag' and 'lead'

### Functions

plotPDF <- function(snap, dat3D, dat2D){
  
  alldata <- data.frame(dat3D,dat2D)
  colnames(alldata)[1] <- '3D'; colnames(lambdas_df)[2] <- '2D'
  #convert data to long format
  lambdas_melt <- melt(alldata, variable.name = "Dimension", value.name = "Lambda")
  
  #set up output plot filename:
  outfn_pdf <- sprintf('snap%04dpdf_%02s.png', snap, proj) # output file name
  #png(filename = file.path(outpath, outfn_pdf), width = 500, height = 500)#, res=40)
  
  #make plot
  plot(
    ggplot(lambdas_melt, aes(x = Lambda, color = Dimension)) +
      geom_line(stat = "density") +
      scale_color_manual(values = c("black", "red")) +
      #lines at median 3D & 2D positions:
      geom_vline(xintercept = median(dat3D), colour = "black") +
      geom_vline(xintercept = median(dat2D), colour = "red", linetype = "dotted") +
      annotate("rect", xmin = median(dat3D) * 0.9, xmax = median(dat3D) * 1.1,
               ymin = -Inf, ymax = Inf, alpha = .2, fill = "blue") + #shaded region 10% around 3D med
      theme_minimal()  +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
      #annotate with p-value & snapshot number:
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = sprintf("p = %4.2e", pval)) +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 3.2, size = 3.6,
               label = ifelse(pval > 1.e-3, TRUE, FALSE)) +
      annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = sprintf("snap %04d", snap))
  )
  #dev.off() #close plot
}

plottimeseries <- function(snap, dat3D, dat2D){
  #Plot lambda time-series showing this snapshot using python script
  
  #concaternate 'quantile' vector to one string:
  quartiles3d <- paste(quantile(lambdas_df$'3D'), collapse = ' ')
  quartiles2d <- paste(quantile(lambdas_df$'2D'), collapse = ' ')
  #concaternate 'quantile' vector to one string:
  sixths3d <- paste(quantile(lambdas_df$'3D', probs = seq(0, 1, 1/6)), collapse = ' ')
  sixths2d <- paste(quantile(lambdas_df$'2D', probs = seq(0, 1, 1/6)), collapse = ' ')
  
  system(paste('python PDFs_lambda.py', outpath, knum, cluster, snap, #model variables
               '--medians3d', medians_3d,
               '--medians2d', medians_2d,
               '--quartiles3d', quartiles3d,
               '--quartiles2d', quartiles2d,
               '--sixths3d', sixths3d,
               '--sixths2d', sixths2d
  ), wait=FALSE)
}

#########
#
# Main:
#

origin <- getwd()

fstr <- 'f16'
qstr <- 'q03'

fn3d <- 'allMSTs_lambar_3D'
proj <- 'xy' #currently doesn't do anything. Change to variable later.
fn2d <- 'allMSTs_lambar_xy'
cluster <- 'cluster_FoV5pc' # use 5 pc field of view

#'root' directory:
rootdir <- '/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0'
knum <- 'k01'

#input (data) and output (plots) directories:
inpath <- file.path(rootdir, paste0(fstr, qstr, '/analysis/runinv_', knum), cluster)
outpath <- file.path(rootdir, 'stats', paste0('pdf_', fstr, qstr, '_', knum))
#create output directory if it doesn't exist:
ifelse(!dir.exists(outpath), dir.create(outpath), FALSE)

nmsts <- 1000 #don't *need* this, it just helps with memory allocation. Might as well, since it's known

nsnaps <- length(readLines(file.path(inpath, 'CDFdata/allMSTs_lambar_3D.dat')))

#Read in lambda data. Row names are from values in first column of table.
all_lambdas3d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_3D.dat'),
                            nrows = nsnaps, row.names = 1, header = FALSE,
                            col.names = c("snap", sprintf("MST%04d", 1:nmsts)))

all_lambdas2d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_xy.dat'),
                            nrows = nsnaps, row.names = 1, header = FALSE,
                            col.names = c("snap", sprintf("MST%04d", 1:nmsts)))

#select snapshots
#snapshots <- c(1:220, 315:349) #k = 1
#snapshots <- c(490:520)  #k = 2
snapshots <- c(300:345)
#snapshots <- c(58, 77, 400, 425)  #k = 3
#snapshots <- c(520:649)  #k = 7
#snapshots <- c(25:42, 60:80, 120:135)  #k = 10

if (nsnaps < tail(snapshots,1)){
  stop("maximum values in 'snapshots' exceeds max snapshot number")
}
compare.spread <- NULL


#get data for each snapshot: nloop lambdas for each snapshot
# Order: [row, column] - data saved for rows specified in 'snapshots', all columns.
#my_lambdas3d <- data.matrix(all_lambdas3d[snapshots, ])
#my_lambdas2d <- data.matrix(all_lambdas2d[snapshots, ])

medians_3d <- data.frame(rowMedians(as.matrix(all_lambdas3d)))
medians_3d <- gsub("c|\\(|\\)|\n|,", "", medians_3d) #remove c, parentheses, \n, and commas from string
medians_2d <- data.frame(rowMedians(as.matrix(all_lambdas2d)))
medians_2d <- gsub("c|\\(|\\)|\n|,", "", medians_2d) #remove c, parentheses, \n, and commas from string


#kernel density estimates
#(see https://www.r-bloggers.com/the-density-function/ for good summary)
for (snapi in 1:length(snapshots)){
  print(snapshots[snapi])
  #'t' function transposes row from table into column for data frame
  lambdas_df <- data.frame(t(all_lambdas3d[snapshots[snapi], ]), t(all_lambdas2d[snapshots[snapi], ]))
  colnames(lambdas_df)[1] <- '3D'; colnames(lambdas_df)[2] <- '2D'
  
  pval <- wilcox.test(as.numeric(lambdas_df$'3D'), as.numeric(lambdas_df$'2D'), paired = FALSE)$p.value
  
  #Call plotting functions defined above:
  #plotPDF(snap = snapshots[snapi], dat3D = lambdas_df$'3D', dat2D = lambdas_df$'2D')
  #plottimeseries(snap = snapshots[snapi], dat3D = lambdas_df$'3D', dat2D = lambdas_df$'2D')
  
  
  #Find y-values at given x-value at various points to find any snaps that don't overlap
  #(to check whether 2D & 3D 'agree' for paper purposes)
  # Get coordinates of density function for 3D & 2D lambdas
  dens3D <- density(lambdas_df$'3D')
  dens2D <- density(lambdas_df$'2D')
  #2D & 3D densities as a function of x (e.g. dens3D_fn(x) returns y for given x):
  dens3D_fn <- approxfun(dens3D$x, dens3D$y)
  dens2D_fn <- approxfun(dens2D$x, dens2D$y)
  
  #For lambda values < than 2 with low p-val ('disagree'), get various density values to compare 2D and 3D.
  if ((median(lambdas_df$'3D') < 2 || median(lambdas_df$'2D') < 2) &&
      (pval < 1.e-3)){
    
    #First, 'NULL' objects in case they were set in the previous iteration:
    # (code works if this isn't done, but can lead to confusion in debugging)
    denspoints_3D <- NULL; denspoints_2D <- NULL
    
    #compare y-values of 2D & 3D density plots at max, min, median, quartiles, & 1/6, 5/6 quantiles:
    denspoints_3D <- c(#dens3D_fn(min(lambdas_df$'3D')),
                       dens3D_fn(quantile(lambdas_df$'3D', .17)),
                       #dens3D_fn(quantile(lambdas_df$'3D', .25)),
                       #dens3D_fn(median(lambdas_df$'3D')),
                       #dens3D_fn(quantile(lambdas_df$'3D', .75)),
                       dens3D_fn(quantile(lambdas_df$'3D', .83))#,
                       #dens3D_fn(max(lambdas_df$'3D'))
    )
    #keep lambdas_df$'3D' for 2D check as x-values need to be the same:
    denspoints_2D <- c(#dens2D_fn(min(lambdas_df$'3D')),
                       dens2D_fn(quantile(lambdas_df$'3D', .17)),
                       #dens2D_fn(quantile(lambdas_df$'3D', .25)),
                       #dens2D_fn(median(lambdas_df$'3D')),
                       #dens2D_fn(quantile(lambdas_df$'3D', .75)),
                       dens2D_fn(quantile(lambdas_df$'3D', .83))#,
                       #dens2D_fn(max(lambdas_df$'3D'))
    )
    
    #replace any 'NA' values with 0:
    denspoints_3D[is.na(denspoints_3D)] <- 0; denspoints_2D[is.na(denspoints_2D)] <- 0
    
    #Find whether the spread of one PDF is encompassed within the other:
    if (all(denspoints_3D >= denspoints_2D) 
        || all(denspoints_3D <= denspoints_2D)){
      
      #NULL any objects that may have been set for a previous snapshot
      x1_3D <- NULL; x2_3D <- NULL; y1 <- NULL; y2 <- NULL
      laglead_df3D <- NULL; laglead_df2D <- NULL
      indices3D <- NULL; indices2D <- NULL
      
      #Either one PDF is within the other or the max tail of one overlaps the min tail of the other.
      #Ensure 3Dx1 > 2Dx1 &  3Dx2 < 2Dx2 or vice versa:
      #3D is easy since the 3D distribution is used as benchmark
      x1_3D <- quantile(lambdas_df$'3D', 0.17)
      x2_3D <- quantile(lambdas_df$'3D', 0.83)
      #2D is trickier.
      #Find the y-values of the relevant quantiles of the 3D distribution:
      y1 <- dens3D_fn(quantile(lambdas_df$'3D', 0.17))
      y2 <- dens3D_fn(quantile(lambdas_df$'3D', 0.83))
      
      #and their corresponding x-values for the 2D distribution;
      #say x1_2D is the first time y1 is reached and x2_2D is the last time y2 is reached.
      
      #First, create a dataframe with | density | density[-1] | density[+1] |
      #using 'lag' and 'lead':
      laglead_df3D <- data.frame(dens = dens3D$y, lag = lag(dens3D$y), lead = lead(dens3D$y))
      laglead_df2D <- data.frame(dens = dens2D$y, lag = lag(dens2D$y), lead = lead(dens2D$y))
      
      #(y1 > lag & y1 < lead) finds indices on positive gradients, (y2 < lag & y2 > lead) for negative
      indices3D <- which(with(laglead_df3D, (y1 > lag & y1 < lead) | (y2 < lag & y2 > lead) ))
      indices2D <- which(with(laglead_df2D, (y1 > lag & y1 < lead) | (y2 < lag & y2 > lead) ))
      
      #If there are more than 4 values (2 either side of y1, y2), stop code
      # & figure out which to take (only implement if needed)
      if (length(indices3D) > 4 || length(indices2D) > 4) {
        stop("2D pdf intersects y1 with positive gradient
             or y2 with negative gradient more than once")
      }
      
      # Take min index as the point to the left of y1 and max index to the right of y2
      indices3D <- c(min(indices3D), max(indices3D))
      indices2D <- c(min(indices2D), max(indices2D))
      #check values:
      indices3D; indices2D
      dens3D$x[indices3D]; dens2D$x[indices2D]
      
      if (
        (dens3D$x[indices3D][1] < dens3D$x[indices2D][1] &&
         dens3D$x[indices3D][2] > dens3D$x[indices2D][2])
        ||
        (dens3D$x[indices3D][1] > dens3D$x[indices2D][1] &&
         dens3D$x[indices3D][2] < dens3D$x[indices2D][2])
      ){
        compare.spread <- c(compare.spread, snapi)
        
        #draw 3D histogram & density function:
        hist(lambdas_df$'3D', freq = F, ylim = c(0,max(dens3D$y,dens2D$y)))
        lines(dens3D, col="red", lwd=2)
        
        abline(v=median(lambdas_df$'3D'), lty=2)
        #median:
        points(median(lambdas_df$'3D'), dens3D_fn(median(lambdas_df$'3D')),
              cex=1.2, pch=20, col="green")
        ##min & max:
        #points(min(lambdas_df$'3D'), dens3D_fn(min(lambdas_df$'3D')),
        #      cex=1.2, pch=20, col="blue")
        #points(max(lambdas_df$'3D'), dens3D_fn(max(lambdas_df$'3D')),
        #      cex=1.2, pch=20, col="blue")
        #1/6 & 5/6 quantiles (0.17, 0.83):
        points(quantile(lambdas_df$'3D', 0.17), dens3D_fn(quantile(lambdas_df$'3D', 0.17)),
              cex=1.2, pch=20, col="green")
        points(quantile(lambdas_df$'3D', 0.83), dens3D_fn(quantile(lambdas_df$'3D', 0.83)),
              cex=1.2, pch=20, col="green")
        ##quartiles (0.17, 0.83):
        #points(quantile(lambdas_df$'3D', 0.25), dens3D_fn(quantile(lambdas_df$'3D', 0.25)),
        #      cex=1.2, pch=20, col="blue")
        #points(quantile(lambdas_df$'3D', 0.75), dens3D_fn(quantile(lambdas_df$'3D', 0.75)),
        #      cex=1.2, pch=20, col="blue")
        
        #draw 2D density function:
        lines(dens2D, col="black", lwd=2)
        text(x=max(dens3D$x,dens2D$x)*0.5,y=max(dens3D$y,dens2D$y),labels=sprintf('snap%04d',snapi))
        #1/6 & 5/6 quantiles (0.17, 0.83):
        points(quantile(lambdas_df$'2D', 0.17), dens3D_fn(quantile(lambdas_df$'3D', 0.17)),
               cex=1.2, pch=4, col="green")
        points(quantile(lambdas_df$'2D', 0.83), dens3D_fn(quantile(lambdas_df$'3D', 0.83)),
               cex=1.2, pch=4, col="green")
        #dev.off()
        
      } #end of 'one within other'
    } #end of 'points of one above points of other' & pval > 0.01
  } #end of 'if median lambdas less than 2'
}
print(compare.spread)

#'ecdf': Empirical Cumulative Distribution Function
#for (snapi in 1:length(snapshots)) {
#  CDF1 <- ecdf(all_lambdas3d[snapshots[snapi], ])
#  CDF2 <- ecdf(all_lambdas2d[snapshots[snapi], ])
#  outfn_cdf <- sprintf('lambdas%02s_snap%04d.png', proj, snapshots[snapi]) # output file name
#  #set up output plot:
#  png(filename = file.path(outpath, outfn_cdf), width = 500, height = 500)#, res=40)
#  plot(CDF1, do.points = FALSE) #plot CDF1
#  lines(CDF2, col='red', do.points = FALSE) #overplot CDF2
#  legend("topleft", legend = c("3D", "xy"), col = c('black', 'red'), lty = c(1, 1))
#  dev.off() #close plot
#}