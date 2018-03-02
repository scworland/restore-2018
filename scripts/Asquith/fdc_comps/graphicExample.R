library(lmomco)
D <- FDClmrdf.nolog
FDC <- D[,c(7:33)]
FF <- names(FDC)
FF <- as.numeric(gsub("f", "", FF))/100
FDC$site <- D$site
FDC$decade <- D$decade
FDC <- FDC[! is.na(FDC$site),]
J <- FDC[FDC$site == "08167000",]
J <- J[J$decade >= 1950,]
qFF <- qnorm(FF);
pdf("junk.pdf", useDingbats=FALSE, width=6.5, height=5)
par(mgp=c(3,0.5,0), las=1)
plot(qFF,J[1,1:27], type="l", log="y", lwd=0.5, col=1,xaxt="n",
     xlab="", ylab="DISCHARGE, CFS", xlim=c(-3.4,3.4), ylim=c(1,50000))
lines(qFF,J[2,1:27], lwd=0.8, col=2)
lines(qFF,J[3,1:27], lwd=1.1, col=4)
lines(qFF,J[4,1:27], lwd=1.4, col=1)
lines(qFF,J[5,1:27], lwd=1.7, col=2)
lines(qFF,J[6,1:27], lwd=2,   col=4)
add.lmomco.axis(las=2, tcl=0.5)
mtext("08167000 Guadalupe River at Comfort, Texas")
dev.off()

