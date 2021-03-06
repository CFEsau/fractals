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
library(forestplot)

#------------
# Functions:
#------------
#================================================================================

plotPDF <- function(snap, dat3D, dat2D){
  
  alldata <- data.frame(dat3D, dat2D)
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
#================================================================================

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
#================================================================================

plotfp <- function(datlabel, medpoint, confint.lo, confint.hi, snapnum, snapcount){
  #Produce a forest plot showing medians & ranges when either
  #the 3D or 2D range falls entirely within the other
  
  #set up output plot filename and forestplot options:
  if(length(snapnum) == 1){
    outfn.fp <- sprintf('ci%02d-%02d_%ssnap%04d.png', #ci for confidence interval
                        as.numeric(sub("%", "", ci[1])),
                        as.numeric(sub("%", "", ci[2])),
                        datlabel[2,1],
                        as.numeric(snapnum)) #use knum & snapshot if just one
    
    plotsize <- c(700, 200)
    xlims <- c(0.5, 2.5)
  } else {
    outfn.fp <- sprintf('ci%02d-%02d_%01d.png',
                        as.numeric(sub("%", "", ci[1])),
                        as.numeric(sub("%", "", ci[2])),
                        snapcount) # use plot number if more than one snapshot
    plotsize <- c(700, 500)
    xlims <- c(0.5, 3.0)
  }
  
  png(filename = file.path(outpath.fp, outfn.fp),
      width = plotsize[1], height = plotsize[2], res = 100)
  
  forestplot(datlabel, medpoint, confint.lo, confint.hi,
             is.summary = c(TRUE, rep(FALSE, length(datlabel)-1)),
             legend = c("3D", "2D"),
             legend_args = fpLegend(pos = list(x = 0.9, y = 0.5),
                                    #legend_args = fpLegend(pos = list(x = 0.95, y = 0.85),
                                    gp = gpar(col = "#CCCCCC", fill = "#F9F9F9")),
             boxsize = 0.1, #size of median point markers
             line.margin = .2,
             col = fpColors(box = c("black", "darkred"), lines = c("black", "darkred")),
             clip = c(0.5, 3.4),
             xticks = c(seq(from = xlims[1], to = xlims[2], by = 0.5)),
             xlab="Lambda",
             txt_gp = fpTxtGp(label = gpar(cex = 0.8),
                              ticks = gpar(cex = 0.8),
                              xlab = gpar(cex = 0.8))
  )
  
  dev.off() #close plot
}
#================================================================================

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

ci <- c("17%", "83%") #confidence interval for forest plots e.g. 17%, 83% for 1/6, 5/6
#ci <- c("25%", "75%")

outpath.fp <- file.path(rootdir, 'stats', 'forestplots')
ifelse(!dir.exists(outpath.fp), dir.create(outpath.fp), FALSE)

#Plot ranges above will vary with data set and quantiles. Automate in future if needed.

for(knum in sprintf("k%02d", 1:10)){
  cat("\nknum: ", knum, "\n") #print knum
  
  #input (data) and output (plots) directories:
  inpath <- file.path(rootdir, paste0(fstr, qstr, '/analysis/runinv_', knum), cluster)
  outpath <- file.path(rootdir, 'stats', paste0('pdf_', fstr, qstr, '_', knum))
  ifelse(!dir.exists(outpath), dir.create(outpath), FALSE)
  
  nmsts <- 1000 #(only used for memory allocation)
  
  nsnaps <- length(readLines(file.path(inpath, 'CDFdata', fn3d)))
  
  #Read in lambda data
  all_lambdas3d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_3D.dat'),
                              nrows = nsnaps, row.names = 1, header = FALSE,
                              col.names = c("snap", sprintf("MST%04d", 1:nmsts)))
  
  all_lambdas2d <- read.table(file.path(inpath, 'CDFdata/allMSTs_lambar_xy.dat'),
                              nrows = nsnaps, row.names = 1, header = FALSE,
                              col.names = c("snap", sprintf("MST%04d", 1:nmsts)))
  
  #select snapshots
  snapshots <- c(1:nsnaps)
  
  if (nsnaps < tail(snapshots, 1)){
    stop("maximum values in 'snapshots' exceeds max snapshot number")
  }
  compare.spread <- NULL #list of snapshots to take a closer look at
  
  medians_3d <- data.frame(rowMedians(as.matrix(all_lambdas3d)))
  medians_3d <- gsub("c|\\(|\\)|\n|,", "", medians_3d) #remove c, parentheses, \n, and commas from string
  medians_2d <- data.frame(rowMedians(as.matrix(all_lambdas2d)))
  medians_2d <- gsub("c|\\(|\\)|\n|,", "", medians_2d) #remove c, parentheses, \n, and commas from string
  
  
  #kernel density estimates
  #(see https://www.r-bloggers.com/the-density-function/ for good summary)
  print("Doing snapshot...")
  for (snapi in 1:length(snapshots)){
    print(snapshots[snapi])
    
    #'t' function transposes row from table into column for data frame
    lambdas_df <- data.frame(t(all_lambdas3d[snapshots[snapi], ]), t(all_lambdas2d[snapshots[snapi], ]))
    colnames(lambdas_df)[1] <- '3D'; colnames(lambdas_df)[2] <- '2D'
    
    pval <- wilcox.test(as.numeric(lambdas_df$'3D'), as.numeric(lambdas_df$'2D'), paired = FALSE)$p.value
    
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
      if((denspoints_3D[ci[1], "x"] < denspoints_2D[ci[1], "x"] && 
          denspoints_3D[ci[2], "x"] > denspoints_2D[ci[2], "x"])
         ||
         (denspoints_3D[ci[1], "x"] > denspoints_2D[ci[1], "x"] && 
          denspoints_3D[ci[2], "x"] < denspoints_2D[ci[2], "x"])
      ){
        
        
        #Call plotting functions defined above:
        #plotPDF(snap = snapshots[snapi], dat3D = lambdas_df$'3D', dat2D = lambdas_df$'2D')
        plottimeseries(snap = snapshots[snapi], dat3D = lambdas_df$'3D', dat2D = lambdas_df$'2D')
        Sys.sleep(1) #pause when plotting time series as python print statements lag behind R
        
        compare.spread <- c(compare.spread, snapshots[snapi])
        
        #Prepare data for forest plot normalised to median 3D value for each PDF
        forestplt[nrow(forestplt)+1, ] <- c(knum, snapshots[snapi], "3D",
                                            denspoints_3D["50%", "x"], #/ denspoints_3D["50%", "x"],
                                            denspoints_3D[ci[1],"x"], #/ denspoints_3D["50%", "x"],
                                            denspoints_3D[ci[2],"x"]) #/ denspoints_3D["50%", "x"])
        forestplt[nrow(forestplt)+1, ] <- c(knum, snapshots[snapi], "2D",
                                            denspoints_2D["50%","x"], #/ denspoints_3D["50%", "x"],
                                            denspoints_2D[ci[1],"x"], #/ denspoints_3D["50%", "x"],
                                            denspoints_2D[ci[2],"x"]) #/ denspoints_3D["50%", "x"])
        
        #draw 3D histogram & density function:
        out.fn <- sprintf('snap%04dpdf_%02s_quantiles.png', snapshots[snapi], proj) # output file name
        #png(filename = file.path(outpath, out.fn), width = 500, height = 500)#, res=40)
        
        #hist(lambdas_df$'3D', freq = F,
        #    xlim = c(min(dens3D$x, dens2D$x), max(dens3D$x, dens2D$x)), ylim = c(0, max(dens3D$y, dens2D$y)))
        plot(dens3D, col="black", lwd=2)
        
        abline(v = median(lambdas_df$'3D'), lty = 2)
        points(denspoints_3D["50%", "x"], denspoints_3D["50%", "y"],
               cex = 1.2, pch = 20, col = "blue") #median
        points(denspoints_3D[ci[1], "x"], denspoints_3D[ci[1], "y"],
               cex = 1.2, pch = 20, col = "blue") # lower quantile / confidence interval
        points(denspoints_3D[ci[2], "x"], denspoints_3D[ci[2], "y"],
               cex = 1.2, pch = 20, col = "blue") # upper quantile / confidence interval
        
        #draw 2D histogram & density function:
        #hist(lambdas_df$'2D', border = "red", freq = F, add=T)
        lines(dens2D, col = "red", lwd = 2)
        text(x = max(dens3D$x, dens2D$x)*0.7, y = max(dens3D$y, dens2D$y),
             labels = sprintf('snap%04d', snapshots[snapi]))
        text(x = max(dens3D$x[1], dens2D$x[1]), y = max(dens3D$y, dens2D$y), labels = knum)
        points(denspoints_2D["50%", "x"], denspoints_2D["50%", "y"],
               cex=1.2, pch=4, col = "blue") #median
        points(denspoints_2D[ci[1], "x"], denspoints_2D[ci[1], "y"],
               cex = 1.2, pch = 4, col = "blue") #lower quantile / confidence interval
        points(denspoints_2D[ci[2], "x"], denspoints_2D[ci[2], "y"],
               cex = 1.2, pch = 4, col = "blue") #upper quantile / confidence interval
        #dev.off()
        
      } #end of 'points of one within points of other'
    } #end of 'if median lambdas less than 2' & pval > 0.01
  } #end of snapshots loop
  cat("List of snapshots for forest plot: ", compare.spread)
} #end of knum loop


forestplt$median <- as.numeric(forestplt$median)
forestplt$lower <- as.numeric(forestplt$lower)
forestplt$upper <- as.numeric(forestplt$upper)


#Set up arguments for 'forestplot' function:
#list of row names; lists of median, lower, and upper values; alignment vector for table columns
ktext <- ifelse(duplicated(forestplt[, 1])[forestplt$dimension == "3D"],
                " ",       #white space if knum is repeated
                forestplt[, 1][forestplt$dimension == "3D"]) #otherwise keep knum label (col. 1)

forestmed <- cbind(forestplt$median[forestplt$dimension == "3D"],
                   forestplt$median[forestplt$dimension == "2D"])
forestlo <- cbind(forestplt$lower[forestplt$dimension == "3D"],
                  forestplt$lower[forestplt$dimension == "2D"])
foresthi <- cbind(forestplt$upper[forestplt$dimension == "3D"],
                  forestplt$upper[forestplt$dimension == "2D"])
forestsnap <- cbind(forestplt$snapshot[forestplt$dimension == "3D"])

if (all(ci == c("17%", "83%"))){
  plotrange <- list(c(1:30), c(31:60), c(61:90), c(91:120), c(121:length(ktext))) #for 1/6, 5/6
} else if (all(ci == c("25%", "75%"))){
  plotrange <- list(c(1:length(ktext)))
} else {
  print("WARNING: no range defined for forest plot. Defaulting to all.")
}

for (i in 1:length(plotrange)){ #'plotrange' is the number of snapshots plotted in
                                # this forest plot to prevent over-crowding
  thisrange <- unlist(plotrange[i])
#  plotfp(datlabel = cbind(c("knum", ktext[thisrange]), c("snap", forestsnap[thisrange])),
#         medpoint = rbind(NA, forestmed[thisrange, ]),
#         confint.lo = rbind(NA, forestlo[thisrange, ]), confint.hi = rbind(NA, foresthi[thisrange, ]),
#         snapnum = forestsnap, snapcount = i)
}
