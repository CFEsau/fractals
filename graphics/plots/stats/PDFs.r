#Produce probability density functions for user-defined snapshots

###########################
#       _ ___      _   _
# |\/| |_  |  |_| / \ | \
# |  | |_  |  | | \_/ |_/
############################

# - Set up directory structure, file name for lambda data, & other variables
# - Read number of lines in fn3d to get number of snapshots
# - Read in lambda data: columns are lambda values, rows are snapshot number
# - Select snapshot range using 'snapshots' variable
# - medians_2d/3d: calculate the median lambda for each row using rowMedians
#       and remove special caracters for passing into plotting function
# Loop over snapshots:
#   - get data frame of 3D & 2D lambda values for this snapshot
#   - perform wilcoxon test to get p value for snapshot
#   - Call plotting functions for PDF & time series
#
#   ########################################################################
#   # Find whether one of the PDFs falls within the extremes of the other:
#   # Pick a couple of points on 3D PDF & see if the x-values of the 2D PDF
#   # for these y values box fall interior/exterior to 3D PDF.
#   - Get (x,y) coordinates at various points on the
#         density functions for the 3D & 2D projections.
#         Saved in dens3D/2D. Get list of coords using e.g. dens3D$x
#   - Get an approximate function for each density curve using approxfun.
#         Saved in dens3D/2D_fn. e.g. dens3D_fn(x) returns y for given x.
# => dens3D/2D are lists of discrete values.
# => dens2D/3D_fn is a continuous interpolated function that takes any value.
#   - Only test for overlap when median lambda of one PDF < 2
#               and p < 1e-3.
#     - Get y-values of various points of interest from 3D plot
#       (e.g. max, min, median, various quantiles)
#     - 



###########################
#           ___
# |\/|  /\   |  |\ |
# |  | /--\ _|_ | \|
############################

library(ggplot2)
library(reshape2)
library(matrixStats) # for rowMedians
library(dplyr) # for 'lag' and 'lead'
 
#------------
# Functions:
#------------
plotPDF <- function(snap, dat3D, dat2D){
  
  alldata <- data.frame(dat3D,dat2D)
  colnames(alldata)[1] <- '3D'; colnames(alldata)[2] <- '2D'
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

#--------
# Main:
#--------

origin <- getwd()

fstr <- 'f16'
qstr <- 'q03'

fn3d <- 'allMSTs_lambar_3D.dat'
proj <- 'xy' #currently doesn't do anything. Change to variable later.
fn2d <- 'allMSTs_lambar_xy.dat'
cluster <- 'cluster_FoV5pc' # use 5 pc field of view

rootdir <- '/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0'

#Set up data frame for error bars that fall within one another:
forestplt <- data.frame(k = character(),
                         snapshot = integer(),
                         dimension = character(),
                         median = double(),
                         lower = double(),
                         upper = double(),
                         stringsAsFactors = FALSE)

for(knum in sprintf("k%02d",1:10)){
  print(knum)

#input (data) and output (plots) directories:
inpath <- file.path(rootdir, paste0(fstr, qstr, '/analysis/runinv_', knum), cluster)
outpath <- file.path(rootdir, 'stats', paste0('pdf_', fstr, qstr, '_', knum))
ifelse(!dir.exists(outpath), dir.create(outpath), FALSE)

nmsts <- 1000 #(only used for memory allocation)

nsnaps <- length(readLines(file.path(inpath, 'CDFdata',fn3d)))

#Read in lambda data
all_lambdas3d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_3D.dat'),
                            nrows = nsnaps, row.names = 1, header = FALSE,
                            col.names = c("snap", sprintf("MST%04d", 1:nmsts)))

all_lambdas2d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_xy.dat'),
                            nrows = nsnaps, row.names = 1, header = FALSE,
                            col.names = c("snap", sprintf("MST%04d", 1:nmsts)))

#select snapshots
snapshots <- c(1:nsnaps)
#snapshots <- c(1:220, 315:349) #k = 1
#snapshots <- c(490:520)  #k = 2
#snapshots <- c(58, 77, 400, 425)  #k = 3
#snapshots <- c(520:649)  #k = 7
#snapshots <- c(25:42, 60:80, 120:135)  #k = 10

if (nsnaps < tail(snapshots,1)){
  stop("maximum values in 'snapshots' exceeds max snapshot number")
}
compare.spread <- NULL #list of snapshots to take a closer look at

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
    
    
    denspoints_3D <- data.frame(x = c(quantile(lambdas_df$'3D', 0), #min
                                      quantile(lambdas_df$'3D', 0.17),
                                      quantile(lambdas_df$'3D', 0.25),
                                      quantile(lambdas_df$'3D', 0.50), #median
                                      quantile(lambdas_df$'3D', 0.75),
                                      quantile(lambdas_df$'3D', 0.83),
                                      quantile(lambdas_df$'3D', 1.0)) #max
    )
    denspoints_3D$y <- dens3D_fn(denspoints_3D$x)
    
    denspoints_2D <- data.frame(x = c(quantile(lambdas_df$'2D', 0), #min
                                      quantile(lambdas_df$'2D', 0.17),
                                      quantile(lambdas_df$'2D', 0.25),
                                      quantile(lambdas_df$'2D', 0.50), #median
                                      quantile(lambdas_df$'2D', 0.75),
                                      quantile(lambdas_df$'2D', 0.83),
                                      quantile(lambdas_df$'2D', 1.0)) #max
    )
    denspoints_2D$y <- dens2D_fn(denspoints_2D$x)
    #If any values are 'NA', stop:
    if(any(is.na(denspoints_3D)) || any(is.na(denspoints_2D))){
      stop("At least one value in denspoints is 'NA'")
    }
    
    
    #Find whether the spread of one PDF is encompassed within the other:
    if((denspoints_3D["17%","x"] < denspoints_2D["17%","x"] && denspoints_3D["83%","x"] > denspoints_2D["83%","x"])
        ||
       (denspoints_3D["17%","x"] > denspoints_2D["17%","x"] && denspoints_3D["83%","x"] < denspoints_2D["83%","x"])
        ){
      
      compare.spread <- c(compare.spread, snapi)
      
      #Prepare data for forest plot normalised to median 3D value for each PDF
      forestplt[nrow(forestplt)+1, ] <- c(knum, snapi, "3D",
                                          denspoints_3D["50%","x"]/denspoints_3D["50%","x"],
                                          denspoints_3D["17%","x"]/denspoints_3D["50%","x"],
                                          denspoints_3D["83%","x"]/denspoints_3D["50%","x"])
      forestplt[nrow(forestplt)+1, ] <- c(knum, snapi, "2D",
                                          denspoints_2D["50%","x"]/denspoints_3D["50%","x"],
                                          denspoints_2D["17%","x"]/denspoints_3D["50%","x"],
                                          denspoints_2D["83%","x"]/denspoints_3D["50%","x"])
        
      #draw 3D histogram & density function:
      #hist(lambdas_df$'3D', freq = F,
       #    xlim = c(min(dens3D$x, dens2D$x), max(dens3D$x, dens2D$x)), ylim = c(0, max(dens3D$y, dens2D$y)))
      plot(dens3D, col="black", lwd=2)
      
      abline(v=median(lambdas_df$'3D'), lty=2)
      points(denspoints_3D["50%","x"], denspoints_3D["50%","y"],
             cex=1.2, pch=20, col="blue") #median
      points(denspoints_3D["17%","x"], denspoints_3D["17%","y"],
             cex=1.2, pch=20, col="blue") #1/6 quantile
      points(denspoints_3D["83%","x"], denspoints_3D["83%","y"],
             cex=1.2, pch=20, col="blue") #5/6 quantile
      
      #draw 2D histogram & density function:
      #hist(lambdas_df$'2D', border = "red", freq = F, add=T)
      lines(dens2D, col="red", lwd=2)
      text(x=max(dens3D$x,dens2D$x)*0.7,y=max(dens3D$y,dens2D$y),labels=sprintf('snap%04d',snapi))
      points(denspoints_2D["50%","x"], denspoints_2D["50%","y"],
             cex=1.2, pch=4, col="blue") #median
      points(denspoints_2D["17%","x"], denspoints_2D["17%","y"],
             cex=1.2, pch=4, col="blue") #median
      points(denspoints_2D["83%","x"], denspoints_2D["83%","y"],
             cex=1.2, pch=4, col="blue") #1/6 quantile
      #dev.off()
      
      # Forest plot to check
      
    } #end of 'points of one within points of other'
  } #end of 'if median lambdas less than 2' & pval > 0.01
} #end of snapshots loop
print(compare.spread)
} #end of knum loop

forestplt$median <- as.numeric(forestplt$median)
forestplt$lower <- as.numeric(forestplt$lower)
forestplt$upper <- as.numeric(forestplt$upper)

fp <- ggplot(data = forestplt[1:40, ], aes(x = as.numeric(median), y=as.numeric(row.names(forestplt[1:40,])),
                                         color = dimension)) +
  scale_color_manual(values = c("black", "red")) +
  geom_errorbarh(aes(xmin = lower, xmax = upper)) +
  geom_point() +
  theme_bw() #+
  #scale_x_continuous(limits=c(0,3), breaks=seq(0,3,by=0.2))
print(fp)


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