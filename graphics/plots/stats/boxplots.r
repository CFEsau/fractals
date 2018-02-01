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

for (f in 1:length(fvals)) {
  for (q in 1:length(qvals)) {
    pdir <- paste0('/local/cfe/backed_up_on_astro3/fractals/r1p0/fbinary0p0/',fstr[f],qstr[q],
                   '/analysis/pvals')
    fn <- file.path(pdir,"pvals.dat")
    if (file.exists(fn)) {
      print(sprintf("f = %.1f, q = %.1f",fvals[f],qvals[q]))
      mydata <- read.table(fn,header=TRUE)
      alldata <- rbind(alldata,mydata)
      } else {
        print(sprintf("%s doesn't exist. Skipping.",file))
        next
    }
  }
}
fbin = "fbinary1p0"
fdim <- c("16", "20", "26", "30")
qvir <- c("03", "05")