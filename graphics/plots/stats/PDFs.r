library(ggplot2)
library(reshape2)

origin <- getwd()

fn3d <- 'allMSTs_lambar_3D'
proj <- 'xy' #currently doesn't do anything. Change to variable later.
fn2d <- 'allMSTs_lambar_xy'
cluster <- 'cluster_FoV5pc'
outpath <- '/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0/f16q03/analysis'

#select snapshots
#knum <- 'k01'
#snapshots <- c(76, 105, 135, 459)
#knum <- 'k02'
#snapshots <- c(200, 275, 400, 420)
#knum <- 'k03'
#snapshots <- c(58, 77, 400, 425)
knum <- 'k10'
snapshots <- c(25:42)#, 60:80, 120:135)

simpath <- file.path(outpath, paste0('runinv_', knum), cluster)
plotsdir <- paste0(outpath, '/plots/pdf_', knum) #outputs directory for plots
#create output directory if it doesn't exist:
ifelse(!dir.exists(plotsdir), dir.create(plotsdir), FALSE)

nmsts <- 1000 #don't *need* this, it just helps with memory allocation. Might as well, since it's known
all_lambdas3d <- read.table(file.path(simpath, 'CDFdata/allMSTs_lambar_3D.dat'), row.names = 1, nrows = nmsts)
all_lambdas2d <- read.table(file.path(simpath, 'CDFdata/allMSTs_lambar_xy.dat'), row.names = 1, nrows = nmsts)

#get data for each snapshot: nloop lambdas for each snapshot
# Order: [row, column] - data saved for rows specified in 'snapshots', all columns.
my_lambdas3d <- data.matrix(all_lambdas3d[snapshots, ])
my_lambdas2d <- data.matrix(all_lambdas2d[snapshots, ])

#kernel density estimates
#(see https://www.r-bloggers.com/the-density-function/ for good summary)
for (snap in 1:length(snapshots)){
  
  df_lambdas <- data.frame(my_lambdas3d[snap, ], my_lambdas2d[snap, ])
  colnames(df_lambdas)[1] <- '3D'; colnames(df_lambdas)[2] <- '2D'
  df_long <- melt(df_lambdas, variable.name = "Dimension", value.name = "Lambda") #convert data to long format
  
  pval <- wilcox.test(as.numeric(my_lambdas3d[snap, ]), as.numeric(my_lambdas2d[snap, ]), paired = FALSE)$p.value
  
  #set up output plot:
  outfn_pdf <- sprintf('pdf%02s_snap%04d.png', proj, snapshots[snap]) # output file name
  png(filename = file.path(plotsdir, outfn_pdf), width = 500, height = 500)#, res=40)
  
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
}


#'ecdf': Empirical Cumulative Distribution Function
#for (snap in 1:length(snapshots)) {
#  CDF1 <- ecdf(my_lambdas3d[snap, ])
#  CDF2 <- ecdf(my_lambdas2d[snap, ])
#  outfn_cdf <- sprintf('lambdas%02s_snap%04d.png', proj, snapshots[snap]) # output file name
#  #set up output plot:
#  png(filename = file.path(plotsdir, outfn_cdf), width = 500, height = 500)#, res=40)
#  plot(CDF1, do.points = FALSE) #plot CDF1
#  lines(CDF2, col='red', do.points = FALSE) #overplot CDF2
#  legend("topleft", legend = c("3D", "xy"), col = c('black', 'red'), lty = c(1, 1))
#  dev.off() #close plot
#}