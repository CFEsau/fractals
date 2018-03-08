
library(ggplot2)
library(Hmisc) # for minor tick marks
library(ggpubr) # for multiplot
library(png) # for reading in & overwriting png at the end

origin <- getwd()

# Data frame for results
agreement_df <- data.frame('fdim'=numeric(),'qvir'=numeric(),'k'=integer(),
                     'U'=numeric(),'t'=numeric(),'gt2'=numeric(),'total'=numeric())

# Build the directory structure:
fdim <- '16'; qvir <- '03' #strings to match directory structure
cluster <- 'cluster_FoV5pc'

# General outputs directory: upper level for all simulation data, data analysis, and plots:
outpath <- paste0('/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0/f',fdim,'q',qvir,'/analysis')
plotsdir <- file.path(outpath,"plots")
pvaldir <- file.path(plotsdir,"pvals")
ifelse(!dir.exists(pvaldir), dir.create(pvaldir), FALSE)

for (k in 1:10) {

  # Directory containing analysed data:
  knum <- ifelse(k<10,paste0('0',as.character(k)),as.character(k))
  simpath <- file.path(outpath,paste0('runinv_k',knum),cluster)
  
  # Input filenames for lambda data:
  fn3d <- 'allMSTs_lambar_3D.dat'; fnxy <- 'allMSTs_lambar_xy.dat'
  # Read in data as data frames:
  df3D <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_3D.dat'), row.names=1)
  dfxy <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_xy.dat'), row.names=1)
  if (ncol(df3D) != ncol(dfxy)) stop("Data frames not of equal size") #check data frame sizes
  nlambdas <- ncol(df3D)  #data frames have 'nloop' columns
  nsnaps <- nrow(df3D)    #and nsnaps rows
  
  # Build list of p-values for each test, and list of median lambdas for 2D & 3D
  pvals <- data.frame();  medianlam <- data.frame()
  
  for (i in 1:nsnaps) {
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
    #med3D <- median(as.numeric(df3D[i,]));  med2D <- median(as.numeric(dfxy[i,]))
    medianlam <- rbind(medianlam,
                       c(
                         median(as.numeric(df3D[i,])), median(as.numeric(dfxy[i,]))
                         )
                       )
    } #end of nsnaps loop
  
  colnames(pvals) <- c( 'U', 't')
  
  #Change p-values of 0 to be really really small
  pvals$U <- ifelse(pvals$U<1.e-99,1.e-99,pvals$U)
  pvals$t <- ifelse(pvals$t<1.e-99,1.e-99,pvals$t)
  
  
  colnames(medianlam) <- c('med_3D', 'med_2D')
  
  #Find whether 3D (& 2D?) median lambdas are above 2.0
  medianlam$agree_median <- ifelse(medianlam$med_3D >= 2.0,# & medianlam$med_2D >= 2.0,
                                #find whether 2D median is within 10% bracket of 3D median
                                ifelse(medianlam$med_2D >= 0.9*medianlam$med_3D & medianlam$med_2D <= 1.1*medianlam$med_3D,
                                       TRUE,
                                       FALSE),
                                #if not both above 2, use p-values
                                'pval')
  medianlam#$pval <- pvals$U
  
  
  #Find p-vales that fall above 0.001 (3 sigma)
  #If true, projection is in agreement with 3D distribution (2D & 3D are same)
  pvals$agree_U <- ifelse(pvals$U>1.e-3,TRUE,FALSE)
  pvals$agree_t <- ifelse(pvals$t>1.e-3,TRUE,FALSE)
  
  #Find projections in agreement using both 10% criterion for lambda>2, and p-value agreement for lambda<2
  medianlam$agree <- ifelse(medianlam$agree_median != 'pval',
                            medianlam$agree_median,
                            pvals$agree_U)
  
  #------------------------------------------------------------
  #
  # Find fraction of lambdas that agree between projections
  #
  nU <- 0; nt <- 0; ngt2 <- 0; ngt2_tot <- 0; ntotal <- 0
  
  #add 1 to nU / nt if agreement is TRUE
  nU <- ifelse(pvals$agree_U,nU+1,nU); nt <- ifelse(pvals$agree_t,nt+1,nt)
  
  #find number of snaphots where lambda3D is above 2.0,
  #and number of times 2D is within 10% of this:
  ngt2 <- ifelse(medianlam$med_2D >= 0.9*medianlam$med_3D & medianlam$med_2D <= 1.1*medianlam$med_3D
                 & medianlam$med_3D >= 2.0, ngt2+1, ngt2)
  
  ngt2_tot <- sum(ifelse(medianlam$med_3D >= 2.0, ngt2_tot+1, ngt2_tot))
  ntotal <- ifelse(medianlam$agree, ntotal+1, ntotal)
  
  #divide number of TRUEs by number of random MSTs to find % agreement 
  pfrac_U <- sum(nU)/nsnaps; pfrac_t <- sum(nt)/nsnaps
  frac_gt2 <- sum(ngt2)/ngt2_tot
  frac_total <- sum(ntotal)/nsnaps
  
  #update knum'th row of data frame:
  agreement_df[as.integer(knum),] <- c(as.numeric(fdim)/10,as.numeric(qvir)/10,
                                    as.integer(knum),pfrac_U,pfrac_t,frac_gt2,frac_total)
  #the above makes 'knum' numeric class. Change column 3 (knum) to 'integer' for formatting
  agreement_df[,3] <- sapply(agreement_df[,3],as.integer)
  
  #==========================
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
agreement_fmt <- agreement_df
agreement_fmt[] <- mapply(sprintf, table_format, agreement_df)

#Output dataframe as table with p-vals to 4 d.p.
write.table(agreement_fmt,file=file.path(pvaldir,"pvals.dat"),
            quote=FALSE, row.names=FALSE, sep="\t")


#box plot of all results across k for this parameter set:
png(filename=paste0(outpath,"/boxplot_f",fdim,"q",qvir,".png"))

boxplot(agreement_df$U,agreement_df$t,agreement_df$gt2,agreement_df$total,
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
