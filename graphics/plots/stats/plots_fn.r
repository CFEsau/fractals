#Used in pvals.r

timeseries_fn <- function(df,xmin=0,xmax=10,dx=1,ymin=0,ymax=15,dy=2,
                          bckcols=c("TRUE" = "green", "FALSE" = "red")){
  

  df_melt <- melt(df, id=c("time_Myr","method","in_agreement"),
                     measure=c("med_3D","med_2D"),
                     variable="dimension", value.name="medianlam") #rename new columns
  
  
  timeseries <- ggplot() +
    geom_line(data=df_melt,aes(x=time_Myr,y=medianlam,group=dimension),
            colour="gray20") +
    labs(x = "Time (Myr)", y = expression(Lambda)) +
    #for colour: geom_line(aes(color=dimension)) +  #colour not relevant - doesn't matter which is 3D/2D
    #scale_color_manual(values=c("brown1", "blue2"))
    geom_rect(data=df,aes(xmin=time_Myr-0.01,xmax=time_Myr,
                                  ymin=-Inf,ymax=ymax,fill=in_agreement),alpha=0.2) +
    scale_color_manual(values=bckcols) +
    theme_minimal() +
    theme(legend.position="none",
          panel.grid.minor = element_blank()) + #remove minor grid lines
    scale_x_continuous(breaks=seq(xmin,xmax,dx)) +
    scale_y_continuous(breaks=seq(ymin,ymax,dy))
  
  #png(filename=paste0(statsdir,"/lambda_",fdim,"q",qvir,"_k",knum,".png"))
  plot(timeseries)
  #dev.off()
}



hist_fn <-function(df,xmax1=1,xmax2=1,xmin1=1e-40,xmin2=1e-12,
                     xlab="log(p)",ylab="nsnaps",ttest=FALSE){
  
  # Wilcoxon/U-test p-val histograms:
  plot1U <- histplot_fn(df,plotdat="U",xmin=xmin1,xmax=xmax1,xlab,ylab)
  plot2U <- histplot_fn(df,plotdat="U",xmin=xmin2,xmax=xmax2,xlab,ylab)
  
  
  if(ttest==TRUE){
    
  # t-test p-val histograms:
  plot1t <- histplot_fn(df,plotdat="t",xmin=xmin1,xmax=xmax1,xlab,ylab)
  plot2t <- histplot_fn(df,plotdat="t",xmin=xmin2,xmax=xmax2,xlab,ylab)
  
  # Add plots (ggarrange adds multiple to same page)
  figure <- ggarrange(plot1U, plot1t, plot2U, plot2t,
                      labels = c("U", "t", "U", "t"),
                      ncol = 2, nrow = 2)
  
  } else {
    
    figure <- ggarrange(plot1U, plot2U, labels = c("U", ""),ncol = 2, nrow = 1)
  }
  #png(filename=paste0(statsdir,"/pvals_",knum,".png"))
  plot(figure)
  #annotate_figure(figure, top=text_grob(paste0("k",knum), face="bold", size=10))
  
  #dev.off()
}

histplot_fn <-function(df,plotdat,xmin,xmax,xlab,ylab){
  
  #countScale.ggplot <- ggplot(p_t$V1, aes(x = p-vals)) +  
  #  geom_histogram() +  
  #  ggtitle("(2) countScale with ggplot") +  
  #  scale_x_log10()
  
  ggplot(data=df,aes_string(x=plotdat)) + geom_histogram(binwidth=1.e0) +
    scale_x_log10(limits = c(xmin, xmax),
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(.x))) + #or format (10^.x)
    #expand=c(0,0)) + #(change close par to comma on above line if uncommenting)
    labs(x = xlab) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() #+
  #annotation_logticks(sides='b')   #log ticks on x-axis
}
