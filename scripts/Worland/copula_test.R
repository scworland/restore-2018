
mu <- 10
sigma <- 2
N <- 1000

x <- rnorm(N,mu,sigma)
pdf <- dnorm(x,mu,sigma)
cdf <- pnorm(x,mu,sigma)

df <- data.frame(x=x,pdf=pdf,cdf=cdf) %>%
  gather(variable,value,-x)

ggplot(df) +
  geom_point(aes(x,value)) +
  facet_wrap(~variable, scales="free",ncol=1) +
  theme_bw()

# copulas
require(mvtnorm)
S <- matrix(c(1,.8,.8,1),2,2) #Correlation matrix
AB <- rmvnorm(mean=c(0,0),sig=S,n=1000) #Our gaussian variables
U <- pnorm(AB) #Now U is uniform - check using hist(U[,1]) or hist(U[,2])
x <- qgamma(U[,1],2) #x is gamma distributed
y <- qbeta(U[,2],1,2) #y is beta distributed
plot(x,y) #They correlate!
