
#--***---***---***---***---***---***---***---***---***---***---***---***---***
#                     DEMONSTRATION OF UNGAGED ESTIMATION
#--***---***---***---***---***---***---***---***---***---***---***---***---***
newdata <- data.frame(x=-50, y=1000, CDA=2.9, ANN_DNI=4.9, MAY=4.8, DEC=4.1)

prob_has_noflow <- plogit(predict(ZERO, newdata=newdata))
pplo <- plogit(predict(PPLO, newdata=newdata))
if(prob_has_noflow < 0.5) pplo <- 0

l1 <- predict(L1, newdata=newdata); l1 <- 10^l1
t2 <- predict(T2, newdata=newdata)
t3 <- predict(T3, newdata=newdata)
t4 <- predict(T4, newdata=newdata)
t5 <- predict(T5, newdata=newdata)
t6 <- predict(T6, newdata=newdata)
lmrvec <- c(l1, l1*t2, t3, t4, t5, t6)
lmr <- vec2lmom(lmrvec)
par <- lmom2par(lmr, type="pe3")

db <- data.frame(x=x, y=y, CDA=CDA, ANN_DNI=ANN_DNI, MAY=MAY, DEC=DEC)
dists <- proxy::dist(db, newdata)
site <- D$site[(1:length(dists))[dists == min(dists)]]
fdc <- D[D$site == site,]
fdc <- as.data.frame(fdc[,7:33])
fdc <- t(fdc); n <- length(fdc); fdc <- fdc[-c(n-1,n),]
ffs <- names(fdc); ffs <- as.numeric(gsub("f", "", ffs))/100

FF <- c(0.0001,seq(0.001,.999,by=.001),0.9999) # convenient nonexceed prob
nFF <- lmomco::f2flo(FF, pp=pplo)
qFF <- qnorm(lmomco::f2f(FF, pp=pplo))
ylim <- range(c(fdc, qlmomco(nFF, par)))
ylim[ylim[1] <= 0] <- c(0.01, ylim[2])
xlim <- qnorm(range(FF))

qdf <- qlmomco(nFF, par)

plot(qFF, qdf, type="l", log="y", xaxt="n",
     xlab="", ylab="QUANTILE, CFS", xlim=xlim, ylim=ylim, col=2, lwd=1.0)
lines(qFF, qdf*L1$duan_smearing, col=4, lwd=1.2)
add.lmomco.axis(las=2, tcl=0.5, side.type="NPP", npp.as.aep=TRUE) # x-axis
mtext("Ungaged FDC estimation (eight regionalization models involved)")
points(qnorm(ffs), fdc)

legend(-4,10000, c("Biased generalized normal",
                   "Unbiased (Duan smearing) generalized normal",
                   "WHA's moonshot (six L-moment fit)"), bty="n", lwd=c(1,1.2,1.4), col=c(2,4,3))
