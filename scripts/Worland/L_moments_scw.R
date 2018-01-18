

# beta 0
b0 <- function(x){
  x <- sort(x)
  n <- length(x)
  temp <- 0
  for(j in 1:n){temp <- temp + x[j]}
  result <- (1/n)*temp
  return(result)
}

# beta 1
b1 <- function(x){
  x <- sort(x)
  n <- length(x)
  temp <- 0
  for(j in 2:n){temp <- temp + (x[j]*((j-1)/(n-1)))}
  result <- (1/n)*temp
  return(result)
}

# beta 2
b2 <- function(x){
  x <- sort(x)
  n <- length(x)
  temp <- 0
  for(j in 3:n){temp <- temp + (x[j]*((j-1)*(j-2))/((n-1)*(n-2)))}
  result <- (1/n)*temp
  return(result)
}

# beta 3
b3 <- function(x){
  x <- sort(x)
  n <- length(x)
  temp <- 0
  for(j in 4:n){temp <- temp + (x[j]*((j-1)*(j-2)*(j-3))/((n-1)*(n-2)*(n-3)))}
  result <- (1/n)*temp
  return(result)
}

# L1
L1 <- function(x){return(b0(x))}

# L2 
L2 <- function(x){return(2*b1(x)-b0(x))}

# L3
L3 <- function(x){return(6*b2(x)-6*b1(x)+b0(x))}

# L4
L4 <- function(x){return(20*b3(x)-30*b2(x)+12*b1(x)-b0(x))}

# t2
T2 <- function(x){return(L2(x)/L1(x))}

# t3
T3 <- function(x){return(L3(x)/L2(x))}

# t4
T4 <- function(x){return(L4(x)/L2(x))}

lmoms_scw <- function(x){return(c(L1(x),L2(x),L3(x),L4(x),T2(x),T3(x),T4(x)))}

# test
x <- rnorm(100,5,3)

obs <- c(lmomco::lmoms(x, nmom=4)$lambdas,lmomco::lmoms(x, nmom=4)$ratios[2:4])
test <- lmoms_scw(x)

obs
test

all.equal(obs,test)


