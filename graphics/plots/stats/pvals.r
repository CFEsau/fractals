library(ggplot2)
library(ggpubr) # for multiplot

origin <- getwd()

# Data frame for results
p_knum <- data.frame('k'=numeric(),'U'=numeric(),'t'=numeric())

# Build the directory structure:
fdim <- '16'; qvir <- '05' #strings to match directory structure
cluster <- 'cluster_FoV5pc'

# General outputs directory: upper level for all simulation data, data analysis, and plots:
outpath <- paste0('/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary1p0/f',fdim,'q',qvir,'/analysis')

for (k in 1:10) {
  # Directory containing analysed data:
  knum <- ifelse(k<10,paste0('0',as.character(k)),as.character(k))
  simpath <- file.path(outpath,paste0('runinv_k',knum),cluster)
  
  plotsdir <- paste0(outpath,'/plots/cdf_',knum) #output directory for plots
  #create output directory if it doesn't exist:
  ifelse(!dir.exists(plotsdir), dir.create(plotsdir), FALSE)
  
  # Input filenames for lambda data:
  fn3d <- 'allMSTs_lambar_3D.dat'; fnxy <- 'allMSTs_lambar_xy.dat'
  # Read in data as data frames:
  df3D <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_3D.dat'), row.names=1)
  dfxy <- read.table(file.path(simpath,'CDFdata/allMSTs_lambar_xy.dat'), row.names=1)
  if (ncol(df3D) != ncol(dfxy)) stop("Data frames not of equal size")
  nlambdas <- ncol(df3D)
  nsnaps <- nrow(df3D)
  
  # Build list of p-values for each test
  pvals <- data.frame() #(p_U=numeric(), p_t=numeric()) didn't work...
  for (i in 1:nsnaps) {
    pvals <- rbind(pvals, c(
      #perform wilcoxon/U-test on data:
      wilcox.test(as.numeric(df3D[i,]), as.numeric(dfxy[i,]))$p.value,
      #perform t-test on data:
      t.test(as.numeric(df3D[i,]), as.numeric(dfxy[i,]))$p.value
      )
    )
  }
  colnames(pvals) <- c( 'U', 't')
  
  #Change p-values of 0 to be really really small
  pvals$U <- ifelse(pvals$U<1.e-99,1.e-99,pvals$U)
  pvals$t <- ifelse(pvals$t<1.e-99,1.e-99,pvals$t)
  
  #Find p-vales that fall below 0.001 (3 sigma)
  #If true, projection is in agreement with 3D distribution
  pvals$agree_U <- ifelse(pvals$U<1.e-3,TRUE,FALSE)
  pvals$agree_t <- ifelse(pvals$t<1.e-3,TRUE,FALSE)
  
  #==========================
  # Wilcoxon/U-test p-val plots:
  
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
  # t-test p-val plots:
  
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
  #plot(figure)
  annotate_figure(figure,
                  top=text_grob("k01",face="bold",size=10)
  )
  #------------------------------------------------------------
  # Find fraction of lambdas that agree between projections
  nU <- 0; nt <- 0
  nU <- ifelse(pvals$agree_U,nU+1,nU); nt <- ifelse(pvals$agree_t,nt+1,nt)
  pfrac_U <- sum(nU)/nlambdas; pfrac_t <- sum(nt)/nsnaps
  
  p_knum[as.numeric(knum),] <- c(as.numeric(knum),pfrac_U,pfrac_t)
  
}

boxplot(p_knum$U,p_knum$t,names=c("U","t"),range=0,ylim=c(0.5,1))