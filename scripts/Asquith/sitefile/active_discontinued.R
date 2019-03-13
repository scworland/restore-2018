library(akqdecay)
#file <- file.choose() # go find full_site_list.csv or decade1950plus_site_list.csv
sites <- read.table(file, header=TRUE, colClasses="character")
sites <- sites$site_no
DV <- new.env()
fill_dvenv(sites, envir=DV, edate="2018-09-30", ignore.provisional=TRUE)
save(DV, file="DVbangP.RData")


beg <- end <- n <- rep(NA, length(ls(DV)))
i <- 0; sites <- sort(ls(DV))
for(i in 1:length(sites)) {
   tmp <- get(sites[i], envir=DV)
   if(length(tmp$year) == 0) {
     message(sites[i]," is missing?")
     next
   }
   yy <- range(tmp$year)
   if(! is.finite(yy[1])) next
   beg[i] <- yy[1]
   end[i] <- yy[2]
   n[i] <- length(tmp$Flow)
}
zz <- data.frame(site_no=sites, count_dvs=n,
                 approx_por_pct=as.integer(100*n/((end-beg)*365.25)),
                 begin_cal_year=beg, end_cal_year=end)
zz$active <- FALSE
zz$active[zz$end_cal_year >= 2017] <- TRUE
zz$active[is.na(zz$count_dvs)] <- FALSE

ExtraSiteFileStuff <- zz
write.table(ExtraSiteFileStuff, file="ExtraSiteFileStuff.txt", quote=TRUE, row.names=FALSE, sep="\t")


