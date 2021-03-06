#To do: Check 'NA' values in tot_agreements. Probably a divide by zero thing.
#       Annotate plots with parameter space/model info

#Calculate p-values
library(data.table) # for melting data frames
library(ggplot2)
library(Hmisc) # for minor tick marks
library(ggpubr) # for multiplot
library(png) # for reading in & overwriting png at the end
source('/local/cfe/backed_up_on_astro3/Github/fractals/graphics/plots/stats/plots_fn.r')


find_agreement <- function(dfcol_method, dfcol_agree, method_type,
                           fdim, qvir, knum, notused){
  
  if(method_type=='all'){
    
    ntot <- length(dfcol_agree) #Count number of snapshots tested, ntot
    # Count number of dfcol_agree where dfcol_agree is TRUE:
    ntot_agree <- length(dfcol_agree[dfcol_agree == 'TRUE'])
    ntot_agree/ntot  #Output fraction of time that projections are in agreement
    
  } else {
    #using a specific method (pval, % tolerance, etc):
    
    #Count the number of times the method is used:
    n_method <- length(dfcol_method[dfcol_method == method_type])
    
    
    #Make a note of the model parameters if this method isn't used
    if (n_method == 0) {
      notused[nrow(notused) + 1, ] <- c(fdim, qvir, knum, method_type)
      #(rbind didn't work... https://stackoverflow.com/questions/5231540/r-losing-column-names-when-adding-rows-to-an-empty-data-frame)
      assign('notused', notused, envir = .GlobalEnv) #equivalent to notused <<- notused (just more explicit)
    }
    
    #Count the number of times the method was *successful*:
    nmethod_agree <- length(dfcol_method[dfcol_method == method_type &
                                           dfcol_agree == 'TRUE'])
    
    #Output fraction of time that lambda for 3D & 2D projection are in agreement:
    nmethod_agree / n_method
    
  }
}

#########
#
#  Main:
#

origin <- getwd()

#Data frame for results: fdim | qvir | k | U | t | in10 | in20 | total
# U & t: % of time p > 0.01 (using all snapshots regardless of lambda value)
# in10: % of time 2D OR 3D are within 10% of each other when one lambda > 2
# in20: % of time 2D AND 3D are within 20% of each other when both lambda > 2
#'total': agreement % for U-test with lambda < 2 and gt10/20 for lambda >= 2
tot_agreement <- data.frame('fdim' = numeric(), 'qvir' = numeric(), 'k' = integer(),
                     'U' = numeric(), 'in10' = numeric(), 'in20' = numeric(), 'total' = numeric())
index <- 0 #row number for tot_agreement

# Variables for directory structure:
fbin <- 'fbinary0p0' ; cluster <- 'cluster_FoV5pc'

fvals <- c(1.6, 2.0, 2.6, 3.0); fstr <- c("f16", "f20", "f26", "f30")
qvals <- c(0.3, 0.5); qstr <- c("q03", "q05")

rootdir <- paste0('/local/cfe/backed_up_on_astro3/fractals/r1p0/', fbin)
#output directory for stats-related plots - create if doesn't exist
statsdir <- file.path(rootdir, 'stats')
ifelse(!dir.exists(statsdir), dir.create(statsdir), FALSE)

pvals_histdir <- file.path(statsdir, 'pvals_hist')
ifelse(!dir.exists(pvals_histdir), dir.create(pvals_histdir), FALSE)

#track instances where lambda doesn't reach 2.0 to cross reference with NaN results in 'agreement':
notused <- data.frame(fdim = double(), qvir = double(), k = integer(),
                      method = character(), stringsAsFactors = FALSE)

for (f in 1:length(fvals)) {
  for (q in 1:length(qvals)) {
    
    modeldir <- paste0(rootdir, '/', fstr[f], qstr[q], '/analysis') #parameter space directory
  
    message(sprintf("\nf = %.1f, q = %.1f", fvals[f], qvals[q])) #print fdim & qvir to console
    
    
    #open plot devices for 'timeseries' and 'method' plots:
    pdf(file = paste0(statsdir, "/lambda_", fstr[f], qstr[q], ".pdf")); timeserdev = dev.cur()
    pdf(file = paste0(statsdir, "/lambda_", fstr[f], qstr[q], "_method.pdf")); methoddev = dev.cur()
    
    for (k in 1:10) {
      #starttime <- Sys.time()
      message(sprintf("\tk = %d", k)) #print knum
      
      knum <- ifelse(k<10, paste0('0', as.character(k)), as.character(k)) #k as string for filepath
      simpath <- file.path(modeldir, paste0('runinv_k', knum), cluster) #data for this model found here
      
      fn3d <- 'allMSTs_lambar_3D.dat'; fnxy <- 'allMSTs_lambar_xy.dat' #Input filenames for lambda data
      
      # Read in data as data frames:
      df3D <- read.table(file.path(simpath, 'CDFdata', fn3d), row.names = 1)
      dfxy <- read.table(file.path(simpath, 'CDFdata', fnxy), row.names = 1)
      if (ncol(df3D) != ncol(dfxy)) stop("Data frames not of equal size") #check data frame sizes
      nlambdas <- ncol(df3D)  #data frames have 'nloop' columns
      nsnaps <- nrow(df3D)    #and nsnaps rows
    
      # Build list of p-values for each test, and list of median lambdas for 2D & 3D
      pvals <- data.frame();  medianlam <- data.frame()
      
      for (i in 1:nsnaps) {
        
        #===========================================
        # Perform U-test & t-test for each snapshot
        #===========================================
        #make list of p-values for all snapshots
        pvals <- rbind(pvals,
                       c( #perform mann-whitney-wilcoxon/U-test on data:
                         wilcox.test(as.numeric(df3D[i, ]), as.numeric(dfxy[i, ]), paired = FALSE)$p.value,
                         #perform t-test on data:
                         t.test(as.numeric(df3D[i, ]), as.numeric(dfxy[i, ]))$p.value
                         )
                       ) #end of pvals rbind
        
        #Get median lambda for 2D & 3D:
        medianlam <- rbind(medianlam,
                           c(median(as.numeric(df3D[i, ])), median(as.numeric(dfxy[i, ])))
                           )
        
        } #end of nsnaps loop
      
      colnames(pvals) <- c( 'U', 't')
      colnames(medianlam) <- c('med_3D', 'med_2D')
      
      #Change p-values of 0 to be really really small
      pvals$U <- ifelse(pvals$U < 1.e-99, 1.e-99, pvals$U)
      pvals$t <- ifelse(pvals$t < 1.e-99, 1.e-99, pvals$t)
      
      #Find p-vales above 0.001 (3 sigma)
      #If true, projection is in agreement with 3D distribution
      pvals$agree_U <- ifelse(pvals$U > 1.e-3, TRUE, FALSE)
      pvals$agree_t <- ifelse(pvals$t > 1.e-3, TRUE, FALSE)
      
      #---------------
      #populate 'medianlam' df with TRUE/FALSE for agreement within 10%/20%:
      
      #make 1D data structures with max median lambda & min median lambda for comparing with 2.0 boundary
      largest <- pmax(medianlam$med_3D, medianlam$med_2D)
      smallest <- pmin(medianlam$med_3D, medianlam$med_2D)
      
      #note on rounding: the 'equals' bit probably won't work due to floating point problems.
      #Instead could use 'round' - e.g. 3*0.8==2.4 gives FALSE, but round(3*0.8,2)==round(2.4,2) gives TRUE.
      #Similarly, could do something with tol <- 1e-2; abs(x-y) <= tol
      #Maybe implement later on, e.g. with #d.p.=3 or 4, or 'tol'=1e-4...
      
      #If the largest of the 2D and 3D values are >= 2.0, see if the smaller lies within 10%:
      medianlam$one_gt2_10 <- ifelse(largest >= 2.0 & smallest < 2.0,
      ifelse(smallest >= 0.9 * largest & smallest <= 1.1 * largest, TRUE, FALSE), 'NA')
      
      #If both 2D and 3D are >= 2.0, see if smaller lies within 20% of larger:
      medianlam$two_gt2_20 <- ifelse(smallest >= 2.0,
      ifelse(smallest >= 0.8 * largest & smallest <= 1.2 * largest, TRUE, FALSE), 'NA')
    
    
      #Combine pvals & medianlam as master data frame with TRUE/FALSE for each snapshot
      # depending on comparison method
      snaps_test <- cbind(pvals, medianlam)
      #if both 2D & 3D median lambda <2, use p-value
      snaps_test$method <- ifelse(largest < 2.0, 'pval',
                                  #if medianlam$one_gt2 is TRUE/FALSE (i.e. not NA), 10% difference is used:
                                  ifelse(medianlam$one_gt2 != 'NA', '10%',
                                       #as above, for 20% difference between lambdas:
                                        ifelse(medianlam$two_gt2 != 'NA', '20%',
                                              'unknown')) #if none of above, 'unknown' (something's wrong!)
                                )
      snaps_test$in_agreement <- ifelse(snaps_test$method == 'pval', snaps_test$agree_U,
                                        ifelse(snaps_test$method == '10%', snaps_test$one_gt2,
                                          ifelse(snaps_test$method == '20%', snaps_test$two_gt2,
                                                 'NA')) #should all be T/F. 'NA' here means something's wrong.
                                 )
      
      # Add 'time_Myr' column to data frame: 0.01 Myr per row
      snaps_test <- cbind(time_Myr = as.numeric(row.names(snaps_test)) * 0.01, snaps_test)
    
      #----------------------------------------------------------#
      # Find fractions of time that projections are in agreement #
      # ----------------------------------------------------------#
      
      #Count the number of times each method is used using 'find_agreement' function.
      #Pass in data frame columns containing method types and T/F, and one method type
      
      frac_p <- find_agreement(snaps_test$method, snaps_test$in_agreement, 'pval',
                               fvals[f], qvals[q], k, notused)
      frac_10pc <- find_agreement(snaps_test$method, snaps_test$in_agreement, '10%',
                                  fvals[f], qvals[q], k, notused)
      frac_20pc <- find_agreement(snaps_test$method, snaps_test$in_agreement, '20%',
                                  fvals[f], qvals[q], k, notused)
      #Total number of snapshots in agreement for this k:
      frac_agree <- find_agreement(snaps_test$method, snaps_test$in_agreement, 'all',
                                   fvals[f], qvals[q], k, notused)
      
      
      #-----------------------
      #Results:
      
      #update next row of data frame with % agreement for each method
      index = index+1
      tot_agreement[index, ] <- c(fvals[f], qvals[q], as.integer(knum),
                                          frac_p, frac_10pc, frac_20pc, frac_agree)
      
      
      ###############################################################
      #                                                             #
      #                        PLOTS                                #
      #                                                             #
      ###############################################################
      
      #===================================
      # Lambda(t) showing 3D & 2D lines with colour coding for T/F:
      
      dev.set(timeserdev) #set plotting device defined earlier
      timeseries_fn(df = snaps_test) #uses median lambdas, saved in both medianlam and snaps_test.
      
      dev.set(methoddev) #set plotting device defined earlier
      methods_fn(df = snaps_test)
      
      #===================================
      # Histograms of p-values from U-test (pval vs. count)
      hist_fn(df = pvals)
      
      # Histograms of p-values from U- and t-tests (pval vs. count)
      #hist_fn(df = pvals, ttest = TRUE)
      
      #endtime <- Sys.time()
      #message(sprintf("\ttime for one k loop = %f s", endtime - starttime)) #time k loop
      } #end of k loop
    
    dev.off(which = timeserdev); dev.off(which = methoddev) #switch off plotting devices
      
      
    #box plot of all results across k for this parameter set
    png(filename = paste0(statsdir, "/boxplot_", fstr[f], qstr[q], ".png"))
      
    #x-axis: U-test box plot, t-test box plot, within 10% & >2 box plot, 'total' in agreement.
    boxplot(tot_agreement$U, tot_agreement$in10, tot_agreement$in20, tot_agreement$total,
            names = c("U", expression(paste(Lambda, "_10%")),
                expression(paste(Lambda, "_20%")),
                    "total"), range = 0, ylim = c(0,1.0))#0.6))
    minor.tick(ny = 4, tick.ratio = 0.3)
    #Grid lines:
    #grid(NA, NULL, col = "lightgray", lty = 2) #nx = NA; ny = NULL (defult major tick positions)
    abline(h = seq(from = 0.0, to = 1.0, by = 0.05), col = "lightgray", lty = 2) #more control than 'grid'
      
    dev.off() #close plot
    
    
    
    ##'annotate_figure' isn't working in the loop - figure not updating.
    ##May just be my version of Rstudio (1.1.383). Could use ggplot instead...
    ##In the meantime, load histograms and annotate the png with knum, then overwrite.
    #for (k in 1:10) {
    #  knum <- ifelse(k < 10, paste0('0', as.character(k)), as.character(k))
    #  
    #  img <- readPNG(paste0(statsdir, "/pvals_", knum, ".png"))
    #  h <- dim(img)[1]
    #  w <- dim(img)[2]
    #  png(paste0(statsdir, "/pvals_", knum, ".png"), width=w, height=h)
    #  par(mar = c(0, 0, 0, 0), xpd = NA, mgp = c(0,0,0), oma = c(0,0,0,0), ann=F)
    #  plot.new()
    #  plot.window(0:1, 0:1)
    #  usr <- par("usr")
    #  rasterImage(img, usr[1], usr[3], usr[2], usr[4] * 0.95)
    #  text(.52, 1.01, paste0("k", knum), cex=1)
    #  
    #  dev.off()
    #}
  
  }#end of qvals loop
}#end of fvals loop


#Format columns for output file.
tot_agreement[, 3] <- sapply(tot_agreement[, 3], as.integer) #Change column 3 (knum) to 'integer'

table_fmt <- c(rep("%.1f", 2), "%2d", rep("%.3f", 4)) #Set formatting

agreement_fmt <- tot_agreement    #Copy df so unformatted one not overwritten
agreement_fmt[] <- mapply(sprintf, table_fmt, tot_agreement) #Make new formatted df


#Output dataframe as table with p-vals to 4 d.p.
#NaN values when there are no values of lambda above the 10% / 20% tolerances (all low lambda)
write.table(agreement_fmt, file = file.path(statsdir, "pvals.dat"),
            quote = FALSE, row.names = FALSE, sep = "\t")

#box plots:
#Read in pvals.dat as tot_agreement if 'allboxplots' needs regenerating)
#melt 'tot_agreement' to fdim | qvir | total
agreement.m <- melt(tot_agreement[, c("fdim", "qvir", "total")], id=c("fdim", "qvir", "total"))

#Create boxplot:
png(filename = paste0(statsdir, "/allboxplots.png"))
plot(
  ggplot(data = agreement.m,                      #plot melted data
         aes(x = as.character(format(fdim, digits = 2)), #x-axis variables formatted to 1 d.p.
             y = total)) +  #y-axis variables
    geom_boxplot(alpha = 0.7,
                 aes(fill = as.character(qvir)), #colour boxes depending on qvir
                 #coef = 20,    #include outliers in whiskers (arbitrary large-ish number)
                 position = position_dodge(width = 0.85), #add space between boxes
                 width = 0.8) + #change width of boxes
    theme_bw() +
    theme(text = element_text(size = 11, family = "Tahoma")) +
    scale_y_continuous(name = "% agreement 2D vs. 3D",
                       breaks = seq(0.0, 1.0, 0.1),    #major ticks
                       limits = c(0.0, 1.0)) +         #y-axis limits
    scale_x_discrete(name = "Fractal dimension") +
    labs(fill = "Virial ratio")
)
dev.off()
