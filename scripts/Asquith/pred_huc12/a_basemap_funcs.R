
setxt1 <- "standard error of fit"
setxt2 <- "Standard error of fit"


east_grids  <- seq(80,100,by=2)
north_grids <- seq(26,38, by=2)

gx <- gy <- vector(mode="numeric")
for(i in 1:length(north_grids)) {
   gy <- c(gy, rep(north_grids[i], length(east_grids)))
   gx <- c(gx, east_grids)
}
GL <- SpatialPoints(cbind(-gx, gy), proj4string=LATLONG)
GL <- spTransform(GL, ALBEA)
XY <- coordinates(GL)
x <- XY[,1]; y <- XY[,2]
#ind <- mgcv::inSide(bnd,x,y)
#XY <- XY[ind,]
GL <- SpatialPointsDataFrame(cbind(-gx, gy), data=data.frame(onoff=rep(1,length(x))),
                                    proj4string=LATLONG)
GL <- spTransform(GL, ALBEA)
ix <- 1:length(x)
#plot(GL, pch=1, col=2)
#text(XY[,1],XY[,2], ix)
GL$onoff[c(1,3:9, 12, 14:20, 23, 34, 45,  56, 67)] <- 0
#rgdal::writeOGR(GL, "gridx2deg", "gridx2deg", driver="ESRI Shapefile")
GL <- GL[-c(1,3:9, 12, 14:20, 23, 34, 45, 56, 67),]
XY <- coordinates(GL)
ix <- 1:length(XY[,1])
#plot(GL, pch=1, col=2)
#text(XY[,1],XY[,2], ix)
GLg <- spTransform(GL, LATLONG)
LL <- coordinates(GLg)
#plot(GL, pch=1, col=2)
#for(i in c(47:56)) {
#  text(XY[i,1], XY[i,2], paste0(as.integer(LL[i,1]),"˚"), cex=0.5, pos=3)
#}
#for(i in c(3,6,16,26,36,46,56)) {
#  text(XY[i,1], XY[i,2], paste0(floor(LL[i,2]+0.001),"˚"), cex=0.5, pos=2)
#}
#for(i in c(1,7,17,27,37,47)) {
#  text(XY[i,1], XY[i,2], paste0(floor(LL[i,2]+0.001),"˚"), cex=0.5, pos=4)
#}

my.north.arrow <-
function(xb, yb, len, lab="North", cex.lab=0.6, tcol="black",  ...) {
    s <- len
    arrow.x = c(-1, 1, 1, 2, 0, -2, -1, -1)
    arrow.y = c(0, 0, 4, 4, 7, 4, 4, 0)
    polygon(xb + arrow.x * s, yb + arrow.y * s, ...)
    text(xb, yb - strheight(lab, cex = cex.lab), lab, cex = cex.lab,
        col = tcol)
}


my.choro.legend <- function(px, py, sh, under="under", over="over", between="to",
                            fmt="%g", cex=1, ...) {
    x = sh$breaks
    lx = length(x)
    if (lx < 3)
        stop("break vector too short")
    res = character(lx + 1)
    res[1] = paste(under, sprintf(fmt, x[1]))
    for (i in 1:(lx - 1)) res[i + 1] <- paste(paste0(sprintf(fmt, x[i]),"1"), # WHA hack
        between, sprintf(fmt, x[i + 1]))
    res[lx + 1] <- paste(over, paste0(sprintf(fmt, x[lx]),"1")) # WHA hack
    maxwidth <- max(strwidth(res))*.75
    temp <- legend(x = px, y = py, legend = rep(" ", length(res)),
        fill = sh$cols, text.width = maxwidth, cex = cex, ...)
    text(temp$rect$left + temp$rect$w, temp$text$y, res, pos = 2,
        cex = cex)
}

my.choro.legend.cat <- function(px, py, sh,
                                fmt="%g", cex=1, cat.levels=NULL, ...) {
    #print(sh)
    x = sh$breaks
    lx = length(x)
    if (lx < 3)
        stop("break vector too short")
    res = cat.levels
    #print(c(length(res), length(sh$cols)))
    maxwidth <- max(strwidth(res))*0.75
    temp <- legend(x = px, y = py, legend = rep(" ", length(res)),
        fill = sh$cols[2:length(res)], text.width = maxwidth, cex = cex, ...)
    offset <- 0 #50000
    text(temp$rect$left + temp$rect$w - offset, temp$text$y, res, pos = 2,
        cex = cex)
}


my.choro.legend.pplo <- function(px, py, sh, under="under", over="over", between="to",
                            fmt="%g", cex=1, ...) {
    under <- "exactly" # WHA hack
    #lx <- length(sh$breaks) # WHA hack
    #sh$breaks <- sh$breaks[2:lx] # WHA hack
    #sh$cols <- sh$cols[2:lx] # WHA hack
    x = sh$breaks
    lx = length(x)
    if (lx < 3)
        stop("break vector too short")
    res = character(lx + 1)
    res[1] = paste(under, sprintf(fmt, x[1]))
    for (i in 1:(lx - 1)) res[i + 1] <- paste(paste0(sprintf(fmt, x[i]),"1"), # WHA hack
        between, sprintf(fmt, x[i + 1]))
    res[lx + 1] <- paste(over, paste0(sprintf(fmt, x[lx]),"1")) # WHA hack
    res[2] <- paste("0 to",x[2]) # WHA hack
    maxwidth <- max(strwidth(res))*.75
    temp <- legend(x = px, y = py, legend = rep(" ", length(res[2:(lx+1)])),
        fill = sh$cols[2:(lx+1)], text.width = maxwidth, cex = cex, ...)
    #print(temp)
    text(temp$rect$left + temp$rect$w, temp$text$y, res[2:(lx+1)], pos = 2,
        cex = cex)
}


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
                "Base modified from USGS digital data, 1:24,000")
  text(145000, 660000, txt, pos=4, cex=0.5)
  plot(GulfStates_modified, add=TRUE, lwd=.4, lty=2)
  STATES <- c("Texas", "Oklahoma", "Missouri", "Arkansas", "Louisiana", "Mississippi",
              "Tennessee", "Kentucky", "Alabama", "Georgia", "Florida")
  STATES <- data.frame(easting=c(-410000, -202900.4,  178000,  178000,  400000,
                                 490000,  740000,  740000,  740000,
                                 1100000, 1290000),
                       northing=c(955139.0, 1400000, 1558716.4, 1400000,  795000,
                                  1165000, 1450000, 1580000, 1325000,
                                  1165000, 800000),
                       state=STATES)
  text(STATES$easting, STATES$northing, STATES$state, pos=4, cex=0.8, col=grey(0.3))
  text(820000, 730000, "Gulf of Mexico", cex=0.9, col=grey(0.3), pos=3)
  plot(GL, lwd=0.4, col=grey(0.22), add=TRUE)
  XY <- coordinates(GL); LL <- coordinates(GLg)
  ladj <- 10000; radj <- 5000; tadj <- 5000
  for(i in c(37:46)) { if(i == 41) next
    txt <- as.integer(LL[i,1])
    text(XY[i,1], XY[i,2]-tadj, paste0(txt,"°"), cex=0.5, pos=3)
    #text(XY[i,1] XY[i,2]-tadj,, paste0(txt,"d"), cex=0.5, pos=3)
  }
  for(i in c(3,6,26,36,46,56)) {
    txt <- floor(LL[i,2]+0.001) # Option-Shift-8 and **NOT** Option-K for the degree
    text(XY[i,1]+ladj, XY[i,2], paste0(txt,"°"), cex=0.5, pos=2)
    #text(XY[i,1]+ladj, XY[i,2], paste0(txt,"d"), cex=0.5, pos=2)
  }
  for(i in c(1,7,17,27,37,47)) {
    txt <- floor(LL[i,2]+0.001)
    text(XY[i,1]-radj, XY[i,2], paste0(txt,"°"), cex=0.5, pos=4)
    #text(XY[i,1]-radj, XY[i,2], paste0(txt,"d"), cex=0.5, pos=4)
  }
  my.north.arrow(1280000, 1380000, 12000, col=grey(0.5))
}

map_base <- function(xlim=NA, ylim=NA) {
  par(lend=1, ljoin=1)
  if(length(grep("spCOV",ls())) != 0) {
    plot(spCOV, pch=NA, xlim=usr[1:2], ylim=usr[3:4])
    plot(GulfStates_modified, add=TRUE, lty=0, col=grey(0.93))
  } else {
    plot(GulfStates_modified, lty=0, col=grey(0.93), xlim=usr[1:2], ylim=usr[3:4])
  }
  polygon(bnd[[1]]$x*1000,bnd[[1]]$y*1000, col=grey(.99), lwd=1.7)
}


map_sebase <- function(data) {
  k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
  for(d in sort(unique(data$decade))) {
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


choropleth_decade_cat <-
  function(data, cuts, n=NA, x=NA) {
    env <- as.environment(as.list(slot(data, "data")))
    x <- get(x, envir=env)
    decade <- get("decade", envir=env)
    cols <- add.alpha(brewer.pal(n+1,"Set2"),.7) # we need one extra color
    k <- 0; ks <- c(0.6,0.8,1.0,1.2,1.4,1.6)
    for(d in sort(unique(decade))) {
      k <- k + 1
      tmp <- x[decade == d]
      shades <- auto.shading(tmp, cutter=cuts, n=n, cols=cols)
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

choropleth_cat <-
  function(data, cuts, n=NA, x=NA, decade="2000") {
    data <- data[data$decade == decade,]
    env <- as.environment(as.list(slot(data, "data")))
    x <- get(x,   envir=env)
    cols <- add.alpha(brewer.pal(n+1,"Set2"),.7)# we need one extra color
    shades <- auto.shading(x, cutter=cuts, n=n+1, cols=cols)
    choropleth(data, x, pch=16, cex=0.4, shading=shades, add=TRUE)
    return(shades)
  }

legend_est <- function(gage="", title="", note=TRUE, shades=NULL,
                       itgage=TRUE, more=NA, sitemap=FALSE,
                       triangle=FALSE, noedwards=TRUE, pplo=FALSE,
                       cat=FALSE, cat.levels=NULL, ...) {
  if(sitemap) {
     note <- FALSE; itgage <- FALSE
  }
  sa  <- paste0("Boundary of study area as defined by Crowley-Ornelas, Knight, and others (2018)")
  tx1 <- paste0("USGS streamgage and symbol* colored according to the legend to the right")
  tx2 <- paste0("Prediction location on NHD+ network colored according to the legend to the right")

  #tx2 <- paste0("'COMID' location in the National Hydrography Dataset\n",
  #              "version 2: Symbol colored by estimated decadal\n",gage)
  #
  xx <- -80000

  if(sitemap) {
    abb <- paste0("Abbreviation: USGS, U.S. Geological Survey")
    text(260000, 290000, abb, cex=0.6, pos=4)
    if(triangle) {
      txt <- paste0("USGS streamgage represented in Crowley-Ornelas, Asquith, Worland, and Knight (2018) ")
    legend(xx, 530000, c(sa, "Edwards aquifer outcrop (Texas Commission on Environmental Quality, 2018)", txt,
        "USGS streamgage with large Edwards aquifer recharge impacts on decadal no-flows and removed from modeling"), bty="n",
           cex=0.7, pt.cex=c(NA,NA,0.8,0.8),
           lwd=c(1.5,3,0.9,0.9), lty=c(1,1,0,0), pch=c(NA,NA,2,2), col=c(1,"#F99B00","#006F41","#8D4200"), ...)
    } else {
      if(noedwards) {
           txt <- paste0("USGS streamgage* represented in Crowley-Ornelas, Asquith, Worland, and Knight (2018) and used in statistical model ")
    legend(xx, 530000, c(sa, "Edwards aquifer outcrop (Texas Commission on Environmental Quality, 2018)", txt), bty="n",
           cex=0.7, pt.cex=c(NA,NA,0.8),
           lwd=c(1.5,3,0.9), lty=c(1,1,0), pch=c(NA,NA,1), col=c(1,"#F99B00","#006F41"), ...)
      } else {
           txt <- paste0("USGS streamgage* represented in Crowley-Ornelas, Asquith, Worland, and Knight (2018) ")
    legend(xx, 530000, c(sa, "Edwards aquifer outcrop (Texas Commission on Environmental Quality, 2018)", txt,
        "USGS streamgage with large Edwards aquifer recharge impacts on decadal no-flows and removed from modeling"), bty="n",
           cex=0.7, pt.cex=c(NA,NA,0.8,0.8),
           lwd=c(1.5,3,0.9,0.9), lty=c(1,1,0,0), pch=c(NA,NA,1,1), col=c(1,"#F99B00","#006F41","#8D4200"), ...)
      }
    text(260000, 350000,
         paste0("* Note that symbol size represents one the six decades:\n     ",
                "1950s (smallest circle) through 2000s (largest circle)"),
         cex=.6, pos=4)
    }
    return()
  }

  if(itgage) {
    legend(xx, 530000, c(sa, tx1, tx2), bty="n", cex=0.7, pt.cex=c(NA,0.8,0.6),
           lwd=c(0.7,0.7,1), lty=c(1,0,0), pch=c(NA,1,16), col=c(1,"#3288BDE6","#D53E4FE6"), ...)
  } else {
    legend(xx, 520000, c(sa, paste0("USGS streamgage and symbol* size changes by decade represented"), tx2), bty="n",
           cex=0.7, pt.cex=c(NA,0.8,0.6),
           lwd=c(0.7,0.7,1), lty=c(1,0,0), pch=c(NA,1,16), col=c(1,grey(0.72),"#D53E4FE6"), ...)
  }
  if(! is.null(shades)) {
   if(pplo) {
      my.choro.legend.pplo(970000, 720000, shades, cex=0.7, bty="n", box.col=grey(1), bg=grey(1),
               fmt="%g", xjust=0, title=title)
   } else if(cat) {
      my.choro.legend.cat(970000, 720000, shades, cex=0.7, bty="n", box.col=grey(1), bg=grey(1),
               fmt="%g", xjust=0, title=title, cat.levels=cat.levels)
   } else {
      my.choro.legend(910000, 720000, shades, cex=0.7, bty="n", box.col=grey(1), bg=grey(1),
               fmt="%g", xjust=0, title=title)
   }
  }
  if(note) {
    text(260000, 380000,
         paste0("* Note that symbol size represents one the six decades:\n     ",
                "1950s (smallest circle) through 2000s (largest circle)"),
         cex=.6, pos=4)
    abb <- paste0("Abbreviations: ",
                  "USGS, U.S. Geological Survey",
                  "; NHD+, National Hydrography Dataset plus version 2,\n     ",
                  "locations by Crowley-Ornelas and others (2018) ")
    if(! is.na(more)) {
      abb <- paste0(abb,more)
      text(260000, 290000, abb, cex=0.6, pos=4)
    } else {
      text(260000, 300000, abb, cex=0.6, pos=4)
    }
  }
}

