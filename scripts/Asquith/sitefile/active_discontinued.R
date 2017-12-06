#sites <- read.table("full_site_list.csv", header=TRUE, colClasses="character")
#sites <- sites$siteno
#DV <- new.env()
#fill_dvenv(sites, envir=DV, edate="2016-09-30", ignore.provisional=TRUE)
#save(DV, file="DVbangP.RData")


beg <- end <- n <- vector(mode="numeric", length(ls(DV)))
i <- 0; sites <- sort(ls(DV))
for(i in 1:length(sites)) {
   tmp <- get(sites[i], envir=DV)
   yy <- range(tmp$year)
   beg[i] <- yy[1]
   end[i] <- yy[2]
   n[i] <- length(tmp$Flow)
}
zz <- data.frame(site_no=sites, count_dvs=n,
                 approx_por_pct=as.integer(100*n/((end-beg)*365.25)),
                 begin_cal_year=beg, end_cal_year=end)
zz$active <- FALSE
zz$active[zz$end_cal_year == 2016] <- TRUE

ExtraSiteFileStuff <- zz
write.table(ExtraSiteFileStuff, file="ExtraSiteFileStuff.txt", quote=TRUE, row.names=FALSE, sep="\t")


