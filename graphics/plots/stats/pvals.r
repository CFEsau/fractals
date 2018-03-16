#Calculate p-values
library(data.table) # for melting data frames
library(ggplot2)
library(Hmisc) # for minor tick marks
library(ggpubr) # for multiplot
library(png) # for reading in & overwriting png at the end

origin <- getwd()

#Set up data frame for results: fdim | qvir | k | U | t | in10 | in20 | total
#U, t, gt2 give agreement frequency: U & t % of time p > 0.01 (using all values),
# in10 gives % of time 2D OR 3D are within 10% of each other when lambda3D > 2.
# in20 gives % of time 2D AND 3D are within 20% of each other when both lambda > 2.
#'total' gives agreement % for U-test with lambda < 2 and gt10/20 for lambda >= 2.
tot_agreement <- data.frame('fdim'=numeric(),'qvir'=numeric(),'k'=integer(),
                     'U'=numeric(), 'in10'=numeric(), 'in20'=numeric(), 'total'=numeric())

# Build the directory structure:
fdim <- '20'; qvir <- '05' #strings to match directory structure
cluster <- 'cluster_FoV5pc'

# General outputs directory: upper level for all simulation data, data analysis, and plots:
outpath <- paste0('/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0/f',fdim,'q',qvir,'/analysis')
#plotsdir <- file.path(outpath,"plots")
pvaldir <- file.path(outpath,"pvals")
ifelse(!dir.exists(pvaldir), dir.create(pvaldir), FALSE)

for (k in 1:10) {

  # Directory containing analysed data:
  knum <- ifelse(k<10,paste0('0',as.character(k)),as.character(k))
  simpath <- file.path(outpath,paste0('runinv_k',knum),cluster)
  
  # Input filenames for lambda data:
  fn3d <- 'allMSTs_lambar_3D.dat'; fnxy <- 'allMSTs_lambar_xy.dat'
  
  # Read in data as data frames:
  df3D <- read.table(file.path(simpath,'CDFdata',fn3d), row.names=1)
  dfxy <- read.table(file.path(simpath,'CDFdata',fnxy), row.names=1)
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
                   c(
                     #perform mann-whitney-wilcoxon/U-test on data:
                     wilcox.test(as.numeric(df3D[i,]), as.numeric(dfxy[i,]),paired=FALSE)$p.value,
                     #perform t-test on data:
                     t.test(as.numeric(df3D[i,]), as.numeric(dfxy[i,]))$p.value
                     )
                   ) #end of pvals rbind
    #Find median lambda for 2D & 3D.
    medianlam <- rbind(medianlam,
                       c(
                         median(as.numeric(df3D[i,])), median(as.numeric(dfxy[i,]))
                         )
                       )
    } #end of nsnaps loop
  
  colnames(pvals) <- c( 'U', 't')
  colnames(medianlam) <- c('med_3D', 'med_2D')
  
  #Change p-values of 0 to be really really small
  pvals$U <- ifelse(pvals$U<1.e-99,1.e-99,pvals$U)
  pvals$t <- ifelse(pvals$t<1.e-99,1.e-99,pvals$t)
  
  #Find p-vales that fall above 0.001 (3 sigma)
  #If true, projection is in agreement with 3D distribution (2D & 3D are same)
  pvals$agree_U <- ifelse(pvals$U>1.e-3,TRUE,FALSE)
  pvals$agree_t <- ifelse(pvals$t>1.e-3,TRUE,FALSE)
  
  #---------------
  
  #save medianlam$med_3D and medianlam$med_2D as 1D data structures to help readability in the next ifelse bit...
  #median3D <- medianlam$med_3D
  #median2D <- medianlam$med_2D
  largest <- pmax(medianlam$med_3D,medianlam$med_2D)
  smallest <- pmin(medianlam$med_3D,medianlam$med_2D)
  
  #note on rounding: the 'equals' bit probably won't work due to floating point problems.
  #Instead could use 'round' - e.g. 3*0.8==2.4 gives FALSE, but round(3*0.8,2)==round(2.4,2) gives TRUE.
  #Similarly, coul do something with tol <- 1e-2; abs(x-y) <= tol
  #Maybe implement later on, e.g. with #d.p.=3 or 4, or 'tol'=1e-4...
  
  #If the largest of the 2D and 3D values are >= 2.0, see if the smaller lies within 10%:
  medianlam$one_gt2_10 <- ifelse(largest >= 2.0 & smallest < 2.0,
    ifelse(smallest >= 0.9*largest & smallest <= 1.1*largest,TRUE,FALSE),'NA')
  
  #If both 2D and 3D are >= 2.0, see if other lies within 20%:
  medianlam$two_gt2_20 <- ifelse(smallest >= 2.0,
    ifelse(smallest >= 0.8*largest & smallest <= 1.2*largest,TRUE,FALSE),'NA')
  
  #test code used in console:
  #a <- c(1.8, 2.1, 2.3, 1.81, 2.3, 2.7, 2.0, 2.4, 2.4)
  #b <- c(1.8, 1.7, 1.9, 2.0, 2.1, 2.4, 1.81, 2.7, 3.0)
  #ifelse((!a >= 2 & b >= 2) | (a >= 2 & !b >= 2),
  #        ifelse((b > 0.9*a & b < 1.1*a) |
  #                   (a > 0.9*b & a < 1.1*b), TRUE, FALSE),
  #        'NA')
  #ifelse(a >= 2 & b >= 2,
  #       ifelse((b > 0.8*a & b < 1.2*a) |
  #                  (a > 0.8*b & a < 1.2*b), TRUE, FALSE),
  #       'NA')
  
  #medianlam <- data.frame(a,b)
  #a_name <- "med_3D"; b_name <- "med_2D"
  #names(medianlam) <- c(a_name,b_name)
  #largest <- pmax(medianlam$med_2D,medianlam$med_3D)
  #smallest <- pmin(medianlam$med_2D,medianlam$med_3D)
  #ifelse((largest >= 2.0 & smallest < 2.0),
  #ifelse((smallest >= 0.9*largest & smallest <= 1.1*largest),TRUE,FALSE),'NA')
  #ifelse(smallest >= 2.0,
  #ifelse(smallest >= 0.8*largest & smallest <= 1.2*largest,TRUE,FALSE),'NA')
  
  #Set up master data frame with TRUE/FALSE for each snapshot depending on comparison method
  snaps_test <- cbind(pvals,medianlam)
  #if both 2D & 3D median lambda <2, use p-value
  snaps_test$method <- ifelse(largest<2.0,'pval',
                              ifelse(medianlam$one_gt2!='NA','10%',
                                     ifelse(medianlam$two_gt2!='NA','20%',
                                            'unknown')) #if none of the above, 'unknown' (something is wrong!)
                              )
  snaps_test$in_agreement <- ifelse(snaps_test$method == 'pval',snaps_test$agree_U,
                                 ifelse(snaps_test$method == '10%',snaps_test$one_gt2,
                                        ifelse(snaps_test$method == '20%',snaps_test$two_gt2,
                                               'NA')) #should all be T/F. An 'NA' here means something is wrong.
                                 )
  
  #----------------------------------------------------------#
  # Find fractions of time that projections are in agreement #
  #----------------------------------------------------------#
  
  #Count the number of times each method is used:
  npval <- length(snaps_test$method[snaps_test$method=='pval'])
  n10pc <- length(snaps_test$method[snaps_test$method=='10%'])
  n20pc <- length(snaps_test$method[snaps_test$method=='20%'])
  #Count the number of times each method was *successful*:
  npval_agree <- length(snaps_test$method[snaps_test$method=='pval' & snaps_test$agree_U=='TRUE'])
  n10pc_agree <- length(snaps_test$method[snaps_test$method=='10%' & snaps_test$one_gt2=='TRUE'])
  n20pc_agree <- length(snaps_test$method[snaps_test$method=='20%' & snaps_test$two_gt2=='TRUE'])
  #Fraction of pvals in within 3sigma (when lambda < 2.0):
  frac_p <- npval_agree/npval
  #Fraction of snapshots within 10/20% of each other:
  frac_10pc <- n10pc_agree/n10pc; frac_20pc <- n20pc_agree/n20pc
  
  #Total number of snapshots in agreement for this k:
  ntot <- length(snaps_test$in_agreement)
  ntot_agree <- length(snaps_test$in_agreement[snaps_test$in_agreement=='TRUE'])
  frac_agree <- ntot_agree/ntot
  
  
  #update knum'th row of data frame:
  tot_agreement[as.integer(knum),] <- c(as.numeric(fdim)/10,as.numeric(qvir)/10,
                                    as.integer(knum),frac_p,frac_10pc,frac_20pc,frac_agree)
  #the above makes 'knum' numeric class. Change column 3 (knum) to 'integer' for formatting
  tot_agreement[,3] <- sapply(tot_agreement[,3],as.integer)
  
  ###############################################################
  #                                                             #
  #                        PLOTS                                #
  #                                                             #
  ###############################################################
  #===================================
  # Lambda(t) showing 3D & 2D lines with colour coding for T/F:
  # Use median lambdas, saved in both medianlam and snaps_test.
  # Add 'time_Myr' column to snaps_test: 0.01 Myr per row
  snaps_test <- cbind(time_Myr=as.numeric(row.names(snaps_test))*0.01,snaps_test)
    
  # Melt 'snapts_test' data frame so 2D are listed under 3D:
  snaps_melt <- melt(snaps_test, id=c("time_Myr","method","in_agreement"),
                      measure=c("med_3D","med_2D"),
                      variable="dimension", value.name="medianlam") #rename new columns
  
  shading <- data.frame(xstart=snaps_test[,"time_Myr"]-0.01, xend=snaps_test[,"time_Myr"],
                        #shade=ifelse(snaps_test[,"in_agreement"]=='TRUE','green',
                        #             ifelse(snaps_test[,"in_agreement"]=='FALSE','red',
                        #                    'gray100')))
                        shade=snapstest[,"in_agreement"])
  
  
  cols <- c("TRUE" = "green", "FALSE" = "red")
  ggplot() +
   geom_line(data=snaps_melt,aes(x=time_Myr,y=medianlam,group=dimension),
             colour="gray20") + 
    #for colour: geom_line(aes(color=dimension)) +  #colour not relevant - doesn't matter which is 3D/2D
    #scale_color_manual(values=c("brown1", "blue2"))
    geom_rect(data=snaps_test,aes(xmin=time_Myr-0.01,xmax=time_Myr,
                               ymin=-Inf,ymax=Inf,fill=in_agreement),alpha=0.2) +
    scale_color_manual(values=cols) +
    theme_minimal() +
    theme(legend.position="none",
          panel.grid.minor = element_blank()) + #remove minor grid lines
    scale_x_continuous(breaks=seq(0,10,1)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=6))
  
  
  
  #===================================
  # Wilcoxon/U-test p-val histograms:
  
  #countScale.ggplot <- ggplot(p_U$V1, aes(x = p-vals)) +  
  #  geom_histogram() +  
  #  ggtitle("(2) countScale with ggplot") +  
  #  scale_x_log10()
  plot1U <- ggplot(data=pvals, aes(pvals$U)) + geom_histogram(binwidth=1.e0) +
    #scale_x_log10(limits = c(1e-12, 1), breaks = seq(min(1e-12),max(1), by = 1.e0))
    scale_x_log10(limits = c(1e-40, 1),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(.x)))+ #or format (10^.x)
                  #expand=c(0,0)) + #(change close par to comma on above line if uncommenting)
    labs(x = "log(p-value), e^-40 -- 1") +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() #+
  #annotation_logticks(sides='b')   #log ticks on x-axis
  
  plot2U <- ggplot(data=pvals, aes(pvals$U)) + geom_histogram(binwidth=1.e0) +
    scale_x_log10(limits = c(1e-12, 1),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(.x))) + #or format (10^.x)
                  #expand=c(0,0)) + #(change close par to comma on above line if uncommenting)
    labs(x = "log(p-value), e^-12 -- 1") +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw()
  
  #==========================
  # t-test p-val histograms:
  
  #countScale.ggplot <- ggplot(p_t$V1, aes(x = p-vals)) +  
  #  geom_histogram() +  
  #  ggtitle("(2) countScale with ggplot") +  
  #  scale_x_log10()
  plot1t <- ggplot(data=pvals, aes(pvals$t)) + 
    geom_histogram(binwidth=1.e0) +
    #scale_x_log10(limits = c(1e-12, 1), breaks = seq(min(1e-12),max(1), by = 1.e0))
    scale_x_log10(limits = c(1e-40, 1),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(.x)))+ #or format (10^.x)
    #expand=c(0,0)) + #(change close par to comma on above line if uncommenting)
    labs(x = "log(p-value), e^-40 -- 1") + 
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() #+
    #annotation_logticks(sides='b')   #log ticks on x-axis
  
  plot2t <- ggplot(data=pvals, aes(pvals$t)) + geom_histogram(binwidth=1.e0) +
    scale_x_log10(limits = c(1e-12, 1),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(.x))) + #or format (10^.x)
    #expand=c(0,0)) + #(change close par to comma on above line if uncommenting)
    labs(x = "log(p-value), e^-12 -- 1") +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw()
  
  # Add plots (ggarrange adds multiple to same page)
  figure <- ggarrange(plot1U, plot1t, plot2U, plot2t,
                 labels = c("U", "t", "U", "t"),
                 ncol = 2, nrow = 2)
  png(filename=paste0(pvaldir,"/pvals_",knum,".png"))
  plot(figure)
  #annotate_figure(figure, top=text_grob(paste0("k",knum), face="bold", size=10))
  
  dev.off()
  
} #end of k loop

#Format columns for output file:
table_format <- c(rep("%.1f",2),"%2d",rep("%.3f",4))
agreement_fmt <- tot_agreement
agreement_fmt[] <- mapply(sprintf, table_format, tot_agreement)

#Output dataframe as table with p-vals to 4 d.p.
write.table(agreement_fmt,file=file.path(pvaldir,"pvals.dat"),
            quote=FALSE, row.names=FALSE, sep="\t")


#box plot of all results across k for this parameter set.
#x-axis: U-test box plot, t-test box plot, within 10% & >2 box plot, 'total' in agreement.
png(filename=paste0(outpath,"/boxplot_f",fdim,"q",qvir,".png"))

boxplot(tot_agreement$U,tot_agreement$t,tot_agreement$gt2,tot_agreement$total,
        names=c("U","t",expression(paste(Lambda,"_2D within 10%")),
                "total"),range=0,ylim=c(0,0.6))

minor.tick(ny=4,tick.ratio=0.3)
#grid(NA,NULL,col="lightgray",lty=2) #nx=NA; ny=NULL (defult major tick positions)
abline(h=seq(from=0.0,to=0.6,by=0.05),col="lightgray",lty=2) #more control than 'grid'

dev.off() #close plot

#'annotate_figure' isn't working in the loop - figure not updating.
#May just be my version of Rstudio (1.1.383). Could use ggplot instead...
#In the meantime, load and annotate the png with knum, then overwrite.
for (k in 1:10) {
  knum <- ifelse(k<10,paste0('0',as.character(k)),as.character(k))
  
  img <- readPNG(paste0(pvaldir,"/pvals_",knum,".png"))
  h <- dim(img)[1]
  w <- dim(img)[2]
  png(paste0(pvaldir,"/pvals_",knum,".png"),width=w,height=h)
  par(mar=c(0,0,0,0),xpd=NA,mgp=c(0,0,0),oma=c(0,0,0,0),ann=F)
  plot.new()
  plot.window(0:1,0:1)
  usr <- par("usr")
  rasterImage(img, usr[1], usr[3], usr[2], usr[4]*0.95)
  text(.52,1.01,paste0("k",knum),cex=1)
  
  dev.off()
}
