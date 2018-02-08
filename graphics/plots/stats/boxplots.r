library(ggplot2)
library(data.table)

#Must run pvals.r first to get data
origin <- getwd()

#Only variable needed is binary fraction.
#Go into fdim & qvir directories, see if pvals.dat exists,
# if it does then read in the data and use to make box plots.
fvals <- c(1.6, 2.0, 2.6, 3.0)
fstr <- c("f16", "f20", "f26", "f30")
qvals <- c(0.3, 0.5)
qstr <- c("q03", "q05")

alldata <- data.table()

fn <- 'alldata.txt'
alldata <- read.table(fn, header=TRUE)
alldata.m <- melt(alldata[,c("qvir","fdim","U")],id=c("qvir","fdim","U"))

#Create boxplot:
plot(ggplot(data=alldata.m,                      #plot melted data
            aes(x=as.character(format(fdim,digits=2)), #x-axis variables formatted to 1 d.p.
                y=U)) +  #y-axis variables
       geom_boxplot(alpha=0.7,
                    aes(fill=as.character(qvir)), #colour boxes depending on qvir
                    coef=20,    #include outliers in whiskers (arbitrary large-ish number)
                    position=position_dodge(width=0.85), #make space between boxes
                    width = 0.8) + #change width of boxes
       theme_bw() +
       theme(text = element_text(size = 11, family = "Tahoma")) +
       scale_y_continuous(name = "% agreement 2D vs. 3D",
                          breaks = seq(0.5, 1.0, 0.1), #major ticks
                          limits=c(0.5,1.0)) +         #y-axis limits
       scale_x_discrete(name = "Fractal dimension") +
       labs(fill = "Virial ratio")
     )
#need plot() around ggplot() for Source... fine without if running from command line...

#box plot using all data combined:
#boxplot(alldata.m$U,range=0,ylim=c(0.5,1))
#minor.tick(ny=4,tick.ratio=0.3)
#abline(h=seq(from=0.5,to=1,by=0.05),col="lightgray",lty=2) #more control than 'grid'

#fbin = "fbinary1p0"
#fdim <- c("16", "20", "26", "30")
#qvir <- c("03", "05")
#write.table(alldata,file="alldata.txt",quote=FALSE,row.names=FALSE,sep="\t")