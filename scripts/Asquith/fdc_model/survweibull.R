library(survival)
load("DEMO.RData")

U <- runif(100) # create of uniform variates
y <- qweibull(U, shape=2, scale=5) # simulate in native R
sm <- survreg(Surv(y)~1, dist="weibull") # fit
yp <- predict(sm) # let R provide the predictions
yp2 <- exp(sm$coefficients[1]) # extract intercept
a <- 1/sm$scale # R's version of the shape factor
b <- exp(sm$coefficients[1]) # R's version of the scale factor

lmr <- lmomco::lmoms(y) # L-moments
print(c(lmr$lambdas[1], lmr$ratios[2]*sqrt(pi)))

x <- qweibull(U, shape=1/a, scale=b)
plot(x,y)
abline(0,1)
plot(y,yp)
abline(0,1)

mu_rweibull <- b * gamma(1 + 1/a)
var_rweibull <- b^2 * (gamma(1 + 2/a) - (gamma(1 + 1/a))^2)
sd_rweibull <- sqrt(var_rweibull)


tZ <-  D[D$decade == "2000",]
tZc <- Surv(log10(tZ$n - tZ$nzero), tZ$nzero != 0, type="right")
tSM <- survreg(tZc~CDA+log10(ppt_mean), data=tZ, dist="weibull")

Yp <- predict(tSM)
R <- residuals(tSM)
Ym <- 0.7115 + 0.0228*tZ$CDA + 0.1729*log10(tZ$ppt_mean)
