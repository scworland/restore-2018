
library(qrnn)

x <- as.matrix(iris[,"Petal.Length",drop=FALSE])
y <- as.matrix(iris[,"Petal.Width",drop=FALSE])

cases <- order(x)
x <- x[cases,,drop=FALSE]
y <- y[cases,,drop=FALSE]

tau <- c(0.05, 0.5, 0.95)

fit1 <- qrnn.fit(x=x, y=y, n.hidden=3, tau=tau[1],iter.max=200, n.trials=1)
qrnn.predict(x,fit1)

set.seed(1)
## QRNN models for conditional 5th, 50th, and 95th percentiles
w <- p <- vector("list", length(tau))
for(i in seq_along(tau)){
  w[[i]] <- qrnn.fit(x=x, y=y, n.hidden=3, tau=tau[i],
                     iter.max=200, n.trials=1)
  p[[i]] <- qrnn.predict(x, w[[i]])
}

## Monotone composite QRNN (MCQRNN) for simultaneous estimation of
## multiple non-crossing quantile functions
x.y.tau <- composite.stack(x, y, tau)
fit.mcqrnn <- qrnn.fit(cbind(x.y.tau$tau, x.y.tau$x), x.y.tau$y,
                       tau=x.y.tau$tau, n.hidden=3, n.trials=1,
                       iter.max=500, monotone=1)
pred.mcqrnn <- matrix(qrnn.predict(cbind(x.y.tau$tau, x.y.tau$x),
                                   fit.mcqrnn), ncol=length(tau))
par(mfrow=c(1, 2))
matplot(x, matrix(unlist(p), nrow=nrow(x), ncol=length(p)), col="red",
        type="l")
points(x, y)
matplot(x, pred.mcqrnn, col="blue", type="l")
points(x, y)

# precip data
data(YVRprecip)
y <- YVRprecip$precip
x <- cbind(sin(2*pi*seq_along(y)/365.25),
           cos(2*pi*seq_along(y)/365.25),
           YVRprecip$ncep)


