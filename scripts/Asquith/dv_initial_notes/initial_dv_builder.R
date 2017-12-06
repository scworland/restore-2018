library(akqdecay) # custom from Asquith built during the project

assign("last.warning", NULL, envir = baseenv()) # purge all warnings

DV <- new.env()
sites <- read.table("starting_siteList_20170202_nowOutdated.txt", header=TRUE, colClasses="character")
sites <- sites$site_no


fill_dvenv(sites, envir=DV, edate="2016-09-30")
save(DV, file="DV.RData")


# See file DVissues.txt, and run those processes.
# Modify DV and then don't forget to save again!

load("DV.RData")


NegativeFlowSites <- c("02300021", "02300082", "02301719", "02310663",
                       "02312700", "02313700", "02322800", "02323592",
                       "02453500", "02462951", "07288955", "07348000",
                       "07353000", "07369000", "07372200", "07380120",
                       "07382500", "08012150", "08041780")

# 1-percent cut off chosen by RRK and WHA for automatic drop from project
for(site in ls(DV)) {
   tmp <- get(site, envir=DV)
   tmp <- tmp[! is.na(tmp$Flow),]
   a <- length(tmp$Flow); b <- length(tmp$Flow[tmp$Flow < 0])
   if(b == 0) next;
   message("#  ", site, " POR=",a, "  NEG=",b, "  PCT=", 100*round(b/a, digits=4))
}
#  02300021 POR=3743  NEG=1263  PCT=33.74
#  02300082 POR=3837  NEG=23  PCT=0.6
#  02301719 POR=4741  NEG=582  PCT=12.28
#  02310663 POR=4493  NEG=826  PCT=18.38
#  02312700 POR=18982  NEG=31  PCT=0.16
#  02313700 POR=15726  NEG=1590  PCT=10.11
#  02322800 POR=5085  NEG=6  PCT=0.12
#  02323592 POR=5512  NEG=3  PCT=0.05
#  02453500 POR=7119  NEG=28  PCT=0.39
#  02462951 POR=13631  NEG=13  PCT=0.1
#  07288955 POR=6574  NEG=35  PCT=0.53
#  07348000 POR=19357  NEG=1  PCT=0.01
#  07353000 POR=8858  NEG=83  PCT=0.94
#  07369000 POR=28308  NEG=133  PCT=0.47
#  07372200 POR=19776  NEG=17  PCT=0.09
#  07380120 POR=8180  NEG=103  PCT=1.26
#  07382500 POR=21644  NEG=62  PCT=0.29
#  08012150 POR=8413  NEG=2078  PCT=24.7
#  08041780 POR=4814  NEG=36  PCT=0.75

RemoveNegSites <- c("02300021", "02301719", "02310663", "02313700",
                    "07380120", "08012150", "08041780")
DV <- akq_rm(RemoveNegSites, envir=DV)
DV <- list2env(DV)

for(site in ls(DV)) {
   tmp <- get(site, envir=DV)
   tmp <- tmp[! is.na(tmp$Flow),]
   a <- length(tmp$Flow); b <- length(tmp$Flow[tmp$Flow < 0])
   if(b == 0) next;
   message("#  ", site, " POR=",a, "  NEG=",b, "  PCT=", 100*round(b/a, digits=4))
   tmp <- get(site, envir=DV)
   tmp$Flow[tmp$Flow < 0] <- NA
   assign(site, tmp, envir=DV)
}
#  02300082 POR=3837  NEG=23  PCT=0.6
#  02312700 POR=18982  NEG=31  PCT=0.16
#  02322800 POR=5085  NEG=6  PCT=0.12
#  02323592 POR=5512  NEG=3  PCT=0.05
#  02453500 POR=7119  NEG=28  PCT=0.39
#  02462951 POR=13631  NEG=13  PCT=0.1
#  07288955 POR=6574  NEG=35  PCT=0.53
#  07348000 POR=19357  NEG=1  PCT=0.01
#  07353000 POR=8858  NEG=83  PCT=0.94
#  07369000 POR=28308  NEG=133  PCT=0.47
#  07372200 POR=19776  NEG=17  PCT=0.09
#  07382500 POR=21644  NEG=62  PCT=0.29

for(site in ls(DV)) {
   tmp <- get(site, envir=DV)
   tmp <- tmp[! is.na(tmp$Flow),]
   a <- length(tmp$Flow); b <- length(tmp$Flow[tmp$Flow < 0])
   if(b == 0) next;
   message("#  ", site, " POR=",a, "  NEG=",b, "  PCT=", 100*round(b/a, digits=4))
}

save(DV, file="DV.RData")
