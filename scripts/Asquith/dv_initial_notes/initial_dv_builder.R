library(akqdecay) # custom from Asquith built during the project

assign("last.warning", NULL, envir = baseenv()) # purge all warnings

DV <- new.env()
sites <- read.table("starting_siteList_20170202_nowOutdated.txt", header=TRUE, colClasses="character")
sites <- sites$site_no


fill_dvenv(sites, envir=DV, edate="2016-09-30", ignore.provisional=FALSE)
save(DV, file="DV.RData")


# See file DVissues.txt, and run those processes.
# Modify DV and then don't forget to save again!

load(file.choose()) # DV.RData


NegativeFlowSites <- c(
 "02300082",
 "02312700",
 "02322800",
 "02323592",
 "02453500",
 "02462951",
 "07288955",
 "07348000",
 "07353000",
 "07369000",
 "07372200",
 "07382500")

# 1-percent cut off chosen by RRK and WHA for automatic drop from project
for(site in ls(DV)) {
   tmp <- get(site, envir=DV)
   tmp <- tmp[! is.na(tmp$Flow),]
   a <- length(tmp$Flow); b <- length(tmp$Flow[tmp$Flow < 0])
   if(b == 0) next;
   message("#  ", site, " POR=",a, "  NEG=",b, "  PCT=", 100*round(b/a, digits=4))
}
#  02300082 POR=3837  NEG=23  PCT=0.6
#  02312700 POR=18982  NEG=31  PCT=0.16
#  02322800 POR=5086  NEG=6  PCT=0.12
#  02323592 POR=5513  NEG=3  PCT=0.05
#  02453500 POR=7119  NEG=28  PCT=0.39
#  02462951 POR=13631  NEG=13  PCT=0.1
#  07288955 POR=6575  NEG=35  PCT=0.53
#  07348000 POR=19357  NEG=1  PCT=0.01
#  07353000 POR=8858  NEG=83  PCT=0.94
#  07369000 POR=28308  NEG=133  PCT=0.47
#  07372200 POR=19776  NEG=17  PCT=0.09
#  07382500 POR=21646  NEG=62  PCT=0.29

RemoveNegSites <- NegativeFlowSites
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


for(site in ls(DV)) {
   tmp <- get(site, envir=DV)
   tmp <- tmp[! is.na(tmp$Flow),]
   b <- tmp$Flow[tmp$Flow < 0]
   if(length(b) == 0) next
   print(site)
   print(tmp$Flow[tmp$Flow < 0])
   tmp$Flow[tmp$Flow < 0] <- 0
   assign(site, tmp, envir=DV)
}

save(DV, file="DVneg2zero.RData")

[1] "02300082"
 [1] -0.45 -0.50 -0.39 -0.23 -0.45 -0.03 -0.80 -0.61 -0.08 -0.36 -0.13 -0.12 -0.60 -0.27
[15] -0.13 -0.04 -0.18 -0.36 -1.24 -0.96 -0.12 -0.16 -0.11
[1] "02312700"
 [1] -20.00 -12.00 -18.00 -19.00 -24.00  -5.14 -22.00 -21.90 -11.70 -56.90 -60.90  -4.00
[13] -28.00 -29.90 -18.00 -11.00  -2.18  -2.81 -17.40  -7.49  -0.22  -0.29  -0.19  -2.01
[25]  -1.96  -0.74  -0.15 -76.00 -78.90 -69.20 -32.40
[1] "02322800"
[1]  -133  -410  -668 -1070 -1020  -558
[1] "02323592"
[1] -335 -935 -197
[1] "02453500"
 [1] -138.0 -151.0 -216.0 -319.0  -32.0  -19.9 -175.0 -221.0 -115.0 -192.0 -172.0 -183.0
[13]  -93.6 -106.0 -121.0 -223.0  -12.1  -19.0  -15.4 -137.0 -354.0 -115.0  -93.5  -35.4
[25] -131.0 -160.0 -111.0  -29.6
[1] "02462951"
 [1]  -5.07  -2.29 -19.90  -2.29  -2.29 -11.10  -4.97  -4.97 -13.40  -4.97 -13.40  -4.97
[13] -30.30
[1] "07288955"
 [1]   -290.0  -1150.0  -2020.0  -2880.0  -3750.0  -2760.0  -1780.0   -790.0   -129.0
[10]    -75.6    -63.7   -336.0  -2920.0  -2920.0  -2310.0  -1790.0  -2760.0  -3480.0
[19]  -4690.0  -7280.0 -14200.0 -18500.0 -24000.0 -32500.0 -41700.0 -46500.0 -49400.0
[28] -48300.0 -33500.0 -36900.0 -37200.0 -40600.0 -37800.0 -24700.0 -10200.0
[1] "07348000"
[1] -50
[1] "07353000"
 [1] -2320.0 -2650.0 -3450.0 -4080.0 -4620.0 -5110.0 -5430.0 -5570.0 -5600.0 -5290.0
[11] -4740.0 -3780.0  -700.0 -1200.0 -2590.0 -3170.0 -3970.0 -3950.0 -3420.0 -2640.0
[21]  -100.0  -600.0 -1110.0 -2670.0 -3780.0 -4700.0 -5390.0 -5980.0 -6370.0 -6570.0
[31] -6200.0 -5180.0 -3550.0  -250.0  -100.0   -42.0  -150.0   -15.0   -25.0   -50.0
[41]   -25.0   -25.0   -25.0   -25.0   -25.0   -25.0   -25.0   -25.0   -25.0   -25.0
[51]   -25.0   -25.0   -50.0   -50.0   -50.0   -50.0   -50.0   -50.0   -50.0   -50.0
[61]   -50.0   -50.0   -50.0   -50.0   -50.0   -50.0   -50.0   -50.0    -5.6  -382.0
[71]  -300.0  -250.0  -150.0  -100.0   -50.0   -50.0   -50.0 -1740.0 -2440.0 -2850.0
[81] -3180.0 -3120.0 -2510.0
[1] "07369000"
  [1] -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200
 [18] -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200
 [35] -200 -200 -200 -100 -100 -100 -100 -100 -200 -200 -200 -200 -200 -200 -200 -200 -200
 [52] -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200 -200
 [69] -200 -200 -200 -200  -50 -150 -400  -50 -150 -200 -350 -200 -100  -50 -100 -300 -400
 [86] -500 -600 -650 -100 -300 -400 -400 -500 -100 -100 -150 -150 -150 -250 -300 -350 -350
[103] -300 -250 -200  -20  -20  -40  -90 -190  -10 -175 -100  -10  -10 -100 -165 -170 -175
[120] -180 -185 -190 -195 -205 -205 -210 -205 -200 -190 -200 -190 -200 -190
[1] "07372200"
 [1]  -4.00 -34.40 -30.00 -26.40 -26.40 -26.40 -25.00 -17.00 -14.00  -8.70  -5.60  -2.91
[13]  -1.60  -1.30  -1.30  -1.11  -0.14
[1] "07382500"
 [1]  -6.00  -3.00 -19.00  -9.40 -11.00 -25.10  -4.51 -19.60  -0.42 -11.70 -48.30 -43.20
[13]  -8.72  -8.13 -28.20  -3.49  -5.36 -43.10 -35.60 -39.90  -8.11 -21.50 -21.80  -7.12
[25] -60.90 -15.90  -7.05 -14.10 -28.80 -21.70  -8.89 -55.40 -15.50  -1.95 -22.90 -13.30
[37]  -8.28 -20.10  -9.38  -7.51 -14.80 -23.90 -15.10 -19.10 -23.20 -16.90 -20.80  -6.56
[49] -15.00 -16.20 -15.80 -15.70  -6.76  -9.76  -5.88  -2.59  -3.82  -4.98 -11.90  -5.44
[61]  -8.51 -18.80
