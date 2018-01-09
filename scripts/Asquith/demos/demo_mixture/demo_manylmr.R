"optimLmoms" <-
function(para, lmr=NULL, minF=0, maxF=1, type1="wei", type2="gev", verbose=FALSE) {
  nmom <- length(para)
  para1 <- lmomco::vec2par(c(para[1], exp(para[2]), para[3]), type=type1)
     if(is.null(para1)) return(Inf)
  para2 <- lmomco::vec2par(c(para[4], exp(para[5]), para[6]), type=type2)
     if(is.null(para2)) return(Inf)
  weight <- 0.5 # pnorm(para[6])
  #print(c(para1$para,para2$para,weight))
  L <- vector(mode="numeric",length=nmom)
  for(r in seq(1,nmom)) { # for each order of moment
    sum <- 0
    for(k in seq(0,r-1)) {
      tmp <- (-1)^k*choose(r-1,k) * exp(lgamma(r+1) - lgamma(r-k-1+1) - lgamma(k+1))
      # Quantile function X(f), which will require numerical integration
      XofF <- function(f) {
            lmomco::par2qua2(f, para1, para2, weight=weight)*f^(r-k-1)*(1-f)^(k) }
      int <- NULL # Perform the numerical integration
      try( int <- integrate(XofF,minF,maxF) ) # if NULL some error detected
      if(is.null(int)) return(Inf)
      sum <- sum + tmp*int$value # sum up
      #if(verbose == TRUE) { # Handy messages
      #  message("abs.err=",int$abs.error, " subdivs=",int$subdivisions,
      #          " text=",int$message)
      #}
    }
    L[r] <- sum/r  # don't forget to divide by the order of the L-moment
  }
  err <- sqrt(sum(sapply(1:nmom, function(i) (L[i] - lmr$lambdas[i])^2)))
  if(verbose) message(" RSE = ",err)
  return(err)
}


library(lmomco)
lmr <- lmomco::vec2lmom(c(305, 263, 0.815, 0.631, 0.473, 0.353))

type1a <- "wei"; type2a <- "gev"
par1 <- lmomco::lmom2par(lmr, type=type1a)
par2 <- lmomco::lmom2par(lmr, type=type2a)
para.init <- c(par1$para[1], log(par1$para[2]), par1$para[3],
               par2$para[1], log(par2$para[2]), par2$para[3])
rt1 <- optim(para.init, fn=optimLmoms, lmr=lmr, type1=type1a, type2=type2a)

type1b <- "gev"; type2b <- "wei"
par1 <- lmomco::lmom2par(lmr, type=type1b)
par2 <- lmomco::lmom2par(lmr, type=type2b)
para.init <- c(par1$para[1], log(par1$para[2]), par1$para[3],
               par2$para[1], log(par2$para[2]), par2$para[3])
rt2 <- optim(para.init, fn=optimLmoms, lmr=lmr, type1=type1b, type2=type2b)


rt <- rt1
par1 <- lmomco::vec2par(c(rt$par[1], exp(rt$par[2]), rt$par[3]), type=type1a)
par2 <- lmomco::vec2par(c(rt$par[4], exp(rt$par[5]), rt$par[6]), type=type2a)
XX1 <- lmomco::par2qua2(FF, para1=par1, para2=par2, weight=0.5)
lmr1 <- lmomco::lmoms(lmomco::par2qua2(runif(100000), para1=par1, para2=par2, weight=0.5))
print(lmr$lambdas)
print(lmr1$lambdas)

rt <- rt2
par1 <- lmomco::vec2par(c(rt$par[1], exp(rt$par[2]), rt$par[3]), type=type1b)
par2 <- lmomco::vec2par(c(rt$par[4], exp(rt$par[5]), rt$par[6]), type=type2b)
XX2 <- lmomco::par2qua2(FF, para1=par1, para2=par2, weight=0.5)
lmr2 <- lmomco::lmoms(lmomco::par2qua2(runif(100000), para1=par1, para2=par2, weight=0.5))
print(lmr$lambdas)
print(lmr2$lambdas)

FF <- lmomco::nonexceeds(); qFF <- qnorm(FF)
ylim <- c(pmax(0.01, min(c(XX1, XX2))), max(c(XX1, XX2)))
plot(qFF, XX1, col=2, type="l", log="y", ylab="QUANTILE", ylim=ylim, lwd=3)
lines(qFF, XX2, col=4)
