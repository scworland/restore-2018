
map_annotation <- function() {
  ss <- list(x=-90000, y=320000)
  SpatialPolygonsRescale(layout.scale.bar(height=.07), offset=ss, scale=200*1000,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  ss <- list(x=-90000, y=303000)
  SpatialPolygonsRescale(layout.scale.bar(height=.08), offset=ss, scale=100*1609.344,
                         fill=c("transparent", "black"), plot.grid=FALSE, lwd=0.4)
  xx <- -122000
  sl <- list(x=xx, y=350000)
  sr <- list(x=xx+200*1000, y=350000)
  text(sl, "0",xx, cex=0.6, pos=4); text(sr, "200 kilometers", cex=0.5, pos=4)
  sl <- list(x=xx, y=283000)
  sr <- list(x=xx+100*1609.344, y=283000)
  text(sl, "0", cex=0.6, pos=4); text(sr, "100 miles", cex=0.5, pos=4)

  txt <- paste0("Albers Equal Area Projection\n",
                "North American Datum of 1983\n",
                "Base modified from U.S. Geological Survey digital data, 1:24,000")
  text(-420000, 1590000, txt, pos=4, cex=0.45)
  plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2)
  STATES <- c("Texas", "Oklahoma", "Missouri", "Arkansas", "Louisiana", "Mississippi",
              "Tennessee", "Kentucky", "Alabama", "Georgia", "Florida")
  STATES <- data.frame(easting=c(-440000, -202900.4,  178000,  178000,  279961.5,
                                 430000,  710000,  710000,  710000,
                                 1100000, 1270000),
                       northing=c(955139.0, 1400000, 1558716.4, 1400000,  795368.5,
                                  1170000, 1450000, 1600000, 1325000,
                                  1170000, 800000),
                       state=STATES)
  text(STATES$easting, STATES$northing, STATES$state, pos=4, cex=0.7, col=grey(0.3))
  plot(GL, lwd=0.4, col=grey(0.22), add=TRUE)
}

map_base <- function(xlim=NA, ylim=NA) {
  par(lend=1, ljoin=1)
  plot(spCOV, pch=NA, xlim=usr[1:2], ylim=usr[3:4])
  plot(GulfStates_modified, add=TRUE, lty=0, col=grey(0.95))
  polygon(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=grey(1), lwd=.7)
}

map_sebase <- function() {
  k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
  for(d in sort(unique(D$decade))) {
    k <- k + 1
    plot(D[D$decade == d,], pch=1, lwd=0.7, cex=ks[k], col=grey(0.72), add=TRUE)
  }
}

choropleth_decade <-
  function(data, cuts, x=NA, rev=FALSE, trans=function(t) {t}) {
    env <- as.environment(as.list(slot(data, "data")))
    if(x == "T2") {
      y <- get("L2", envir=env)
      x <- get("L1", envir=env)
      x <- y/x
    } else {
      x <- get(x, envir=env)
    }
    decade <- get("decade", envir=env)
    cols <- add.alpha(brewer.pal(10,"Spectral"),.7)
    if(rev) cols <- rev(cols)
    k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
    for(d in sort(unique(decade))) {
      k <- k + 1
      tmp <- trans(x[decade == d])
      shades <- auto.shading(tmp, cutter=cuts, n=9, cols=cols)
      choropleth(data[data$decade == d,], tmp, pch=1, lwd=0.7, cex=ks[k], shading=shades, add=TRUE)
    }
  }

choropleth_cov <-
  function(data, cuts, x=NA, decade="2000", rev=FALSE, trans=function(t) {t}) {
    data <- data[data$decade == decade,]
    env <- as.environment(as.list(slot(data, "data")))
    x <- get(x,   envir=env)
    cols <- add.alpha(brewer.pal(10,"Spectral"),.7)
    if(rev) cols <- rev(cols)
    tmp <- trans(x)
    shades <- auto.shading(tmp, cutter=cuts, n=9, cols=cols)
    choropleth(data, tmp, pch=16, cex=0.4, shading=shades, add=TRUE)
    return(shades)
  }

legend_est <- function(gage="", title="", note=TRUE, shades=NA, itgage=TRUE, ...) {
  tx1 <- paste0("U.S. Geological Survey streamflow-gaging station:\n",
                "Symbol colored* by observed decadal ",gage)
  tx2 <- paste0("'COMID' location in the National Hydrography Dataset\n",
                "version 2: Symbol colored by estimated decadal\n",gage)
  if(itgage) {
    legend(-80000, 595000, c(tx1, tx2), bty="n", cex=0.7, pt.cex=c(0.8,0.6),
           lwd=c(0.7,1), lty=c(0,0), pch=c(1,16), col=c("#3288BDE6","#D53E4FE6"), ...)
  } else {
    legend(-80000, 520000, c(tx2), bty="n", cex=0.7, pt.cex=c(0.6),
           lwd=c(1), lty=c(0), pch=c(16), col=c("#D53E4FE6"), ...)
  }
  choro.legend(850000, 690000, shades, cex=0.7, bty="n", box.col=grey(1), bg=grey(1),
               fmt="%g", xjust=0, title=title)
  if(note) {
    text(270000, 320000,
         "* Note that symbol size represents decade\n (1950 [smallest circle] through 2000 [largest circle]).",
         cex=.6, pos=4)
  }
}

setxt1 <- "standard error of fit"
setxt2 <- "Standard error of fit"

