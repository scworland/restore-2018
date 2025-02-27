# source("a_pred_huc12.R")
source("a_basemap_funcs.R")


#-----------------------------------------------------------------------
pdf("PPLOfit_junk.pdf", useDingbats=FALSE, width=11, height=10)
  plot(spRESTORE_MGCV_BND)  # by creation of the PDF, we can get a handle on a global
  usr <- par()$usr # setting of the plotting limits by preserving the usr.
dev.off()
unlink("PPLOfit_junk.pdf")  # just quietly throw the file away

manypdfs <- TRUE
file <- "PPLOfit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="pplo", cuts=pploCuts, rev=TRUE)
    shades <- choropleth_cov(H12PPLOdf, decade=d, x="est_pplo", cuts=pploCuts, rev=TRUE)
    legend_est(gage="no-flow fraction",
               title=paste0(d," decade\n","no-flow fraction"),
               note=TRUE, shades=shades, pplo=TRUE)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
file <- "PPLOsefit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase(D)
    shades <- choropleth_cov(H12PPLOdf, decade=d, x="se.fit_flowtime", cuts=pploCutsSE, rev=TRUE)
    legend_est(gage=setxt1,
               title=paste0(d," decade\n",setxt1),
               note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
#-----------------------------------------------------------------------
file <- "L1fit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  H12L1df$est_L1_log10 <- log10(H12L1df$est_L1)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="L1", cuts=L1Cuts, trans=log10)
    shades <- choropleth_cov(H12L1df, decade=d, x="est_L1_log10", cuts=L1Cuts)
    legend_est(gage="mean nonzero streamflow",
               title=paste0(d," decade\n","mean nonzero streamflow,\n","in log10(cms)"),
               note=TRUE, shades=shades, more="; cms, cubic meters per second")
    map_annotation()
    if(manypdfs) dev.off()
  }
  H12L1df$est_L1_log10 <- NULL
if(! manypdfs) dev.off()
file <- "L1sefit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase(D)
    shades <- choropleth_cov(H12L1df, decade=d, x="se.fit_L1", cuts=L1CutsSE, rev=TRUE)
    legend_est(gage=setxt1,
               title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
#-----------------------------------------------------------------------
file <- "T2fit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4]);
    choropleth_decade(D, x="T2", cuts=T2Cuts)
    shades <- choropleth_cov(H12T2df, decade=d, x="est_T2", cuts=T2Cuts)
    legend_est(gage="LCV of nonzero streamflow",
               title=paste0(d," decade\n","L-CV of nonzero streamflow"),
               note=TRUE, shades=shades, more="; L-CV, coefficient of L-variation")
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
file <- "T2sefit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase(D)
    shades <- choropleth_cov(H12T2df, decade=d, x="se.fit_T2", cuts=T2CutsSE, rev=TRUE)
    legend_est(gage=setxt1,
               title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
#-----------------------------------------------------------------------
file <- "T3fit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="T3", cuts=T3Cuts)
    shades <- choropleth_cov(H12T3df, decade=d, x="est_T3", cuts=T3Cuts)
    legend_est(gage="L-skew of nonzero streamflow",
               title=paste0(d," decade\n","L-skew of nonzero streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
    #unique(spCOV$comid[abs(spCOV$dec_long_va - -97.92399) < .1 &
    #                   abs(spCOV$dec_lat_va  -  29.56245) < .1])
    xy <- coordinates(spCOV[spCOV$comid == "1620855", ])
    points(xy[,1], xy[,2], pch=0, col="#b51b96", cex=2)
    arrows(x0=147500, y0=573000, xy[,1], xy[,2],
         lwd=.5, length=0.15, angle=15, col="#b51b96")
    text(147500, 573000, "COMID 1620855 used for a discussion point in the text",
                       cex=0.7, pos=4, col="#7a1265")
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
file <- "T3sefit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase(D)
    shades <- choropleth_cov(H12T3df, decade=d, x="se.fit_T3", cuts=T3CutsSE, rev=TRUE)
    legend_est(gage=setxt1,
               title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
#-----------------------------------------------------------------------
file <- "T4fit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="T4", cuts=T4Cuts)
    shades <- choropleth_cov(H12T4df, decade=d, x="est_T4", cuts=T4Cuts)
    legend_est(gage="L-kurtosis of nonzero streamflow",
               title=paste0(d," decade\n","L-kurtosis of nonzero streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
file <- "T4sefit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase(D)
    shades <- choropleth_cov(H12T4df, decade=d, x="se.fit_T4", cuts=T4CutsSE, rev=TRUE)
    legend_est(gage=setxt1,
               title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
#-----------------------------------------------------------------------
file <- "T5fit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="T5", cuts=T5Cuts)
    shades <- choropleth_cov(H12T5df, decade=d, x="est_T5", cuts=T5Cuts)
    legend_est(gage="Tau5 of nonzero streamflow",
               title=paste0(d," decade\n","Tau5 of nonzero streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()
file <- "T5sefit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4]); map_sebase(D)
    shades <- choropleth_cov(H12T5df, decade=d, x="se.fit_T5", cuts=T5CutsSE, rev=TRUE)
    legend_est(gage=setxt1,
               title=paste0(d," decade\n",setxt1), note=FALSE, shades=shades, itgage=FALSE)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()


D$overL1 <- (1-D$pplo)*(D$L1)
file <- "OverL1fit"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  OverL1df$est_overL1_log10 <- log10(OverL1df$est_overL1)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="overL1", cuts=OverL1Cuts, trans=log10)
    shades <- choropleth_cov(OverL1df, decade=d, x="est_overL1_log10", cuts=OverL1Cuts)
    legend_est(gage="mean streamflow",
               title=paste0(d," decade\n","mean streamflow,\n","in log10(cms)"),
               note=TRUE, shades=shades, more="; cms, cubic meters per second")
    map_annotation()
    if(manypdfs) dev.off()
  }
  OverL1df$est_overL1_log10 <- NULL
if(! manypdfs) dev.off()
D$overL1 <- NULL



spCOV$special_spatial <- retransin(spCOV$grassland)/100
D$special_spatial <- retransin(D$grassland)/100
file <- "GrassLand"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="special_spatial", cuts=grassCuts)
    shades <- choropleth_cov(spCOV, decade=d, x="special_spatial", cuts=grassCuts)
    legend_est(gage="fraction grassland",
               title=paste0(d," decade\n","fraction grassland"),
               note=TRUE, shades=shades)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()


spCOV$special_spatial <- retransin(spCOV$developed)/100
D$special_spatial <- retransin(D$developed)/100
file <- "Developed"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="special_spatial", cuts=developedCuts)
    shades <- choropleth_cov(spCOV, decade=d, x="special_spatial", cuts=developedCuts)
    legend_est(gage="fraction developed",
               title=paste0(d," decade\n","fraction developed"),
               note=TRUE, shades=shades)
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()



bedpermCuts <- function(x, n=6, ...) {
   labs <- 1:n
   cuts <- labs
   cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

spCOV$special_spatial <- as.numeric(spCOV$bedperm)
D$special_spatial <- as.numeric(D$bedperm)
file <- "BedPerm"
if(! manypdfs) pdf(paste0(file,".pdf"), useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(manypdfs) pdf(paste0(file,"_",d,".pdf"), useDingbats=FALSE, width=11, height=10)
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade_cat(D, x="special_spatial", cuts=bedpermCuts, n=6)
    shades <- choropleth_cat(spCOV, decade=d, x="special_spatial", cuts=bedpermCuts, n=6)
    legend_est(gage="bed permeability class",
               title=paste0(d," decade\n","bed permeability class"),
               note=TRUE, shades=shades, cat=TRUE, cat.levels=levels(D$bedperm))
    map_annotation()
    if(manypdfs) dev.off()
  }
if(! manypdfs) dev.off()



#------------------------------------------------------------------
#
#
#
#
quantile(H12L1df$delta_est_L1, probs=(1:9)/10, na.rm=TRUE)

L1delCuts <- function(x, n=9, ...) {
  labs <- 1:n
  cuts <- c(-1, -.5, -.2, -0.05, 0, 0.05, 0.2, 0.5, 1)
  cuts <- cuts[labs]; names(cuts) <- paste("#", labs, sep=""); cuts
}

pdf("L1del.pdf", useDingbats=FALSE, width=11, height=10)
  for(d in sort(unique(D$decade))) {
    if(d == "1950") next
    map_base(xlim=usr[1:2], ylim=usr[3:4])
    choropleth_decade(D, x="L1", cuts=L1delCuts)
    shades <- choropleth_cov(H12L1df, decade=d, x="delta_est_L1", cuts=L1delCuts)
    legend_est(gage="L1 change of nonzero streamflow", title=paste0(d," decade\n","Change in L1 of nonzero streamflow"),
               note=TRUE, shades=shades)
    map_annotation()
  }
dev.off()

cropthem <- TRUE
crop_em <- function(decades=NULL, spawn=FALSE, files=NULL) {
  if(is.null(files)) {
     files <- c("PPLOfit", "L1fit", "OverL1fit", "T2fit", "T3fit", "T4fit", "T5fit")
     files <- c(files, "GrassLand", "BedPerm", "Developed")
     files <- c(files, "PPLOsefit", "L1sefit", "T2sefit", "T3sefit",
                                               "T4sefit", "T5sefit")
  }
  for(file in files) { for(d in decades) { my.file <- paste0(file,"_",d,".pdf")
    system(paste0("pdfcrop --margins '-46 -110 -43 0' --clip ",my.file," ",my.file))
  }}
}
crop_em(decades=sort(unique(D$decade)), spawn=cropthem)



