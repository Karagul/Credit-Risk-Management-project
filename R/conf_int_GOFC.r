library(MKmisc)
library(binom)


a <- 0.2 ## correlation between assets (global var)

###### gofc_prob calculates the probability of d defaults in gofc

gofc_prob <- function(d,n,p) {
f <- function(m) {pnorm((qnorm(p)-a*m)/sqrt(1-a^2))^d*(1-pnorm((qnorm(p)-a*m)/sqrt(1-a^2)))^(n-d)*dnorm(m)}
return(choose(n,d) * integrate(f,-Inf, Inf)$value)
}

###### vectorize gofc_prob

gofc_prob_vec <- function(n,p) {
x <- rep(0,n+1)
for (i in 0:n) x[i+1] <- gofc_prob(i,n,p)
return(x)
}

clopper_gofc <- function(n,p,conf.level=.95){
L <- U <- rep(0,n+1)
x <- 0:n   
pt1 <- qnorm(qbeta((1-conf.level)/2, x, n-x+1))*sqrt(1-a^2)
pt2 <- qnorm(qbeta(1-(1-conf.level)/2, x+1, n-x))*sqrt(1-a^2)
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt1[i] + a*m)*dnorm(m)}
g <- function(m) {pnorm(pt2[i] + a*m)*dnorm(m)}
L[i] <- integrate(f,-Inf,Inf)$value   
U[i] <- integrate(g,-Inf,Inf)$value  
}
L[1]<-0
U[n+1]<-1
return(c(L,U))
#return(c(cp(n,p,L,U),el(n,p,L,U)))
}


wald_gofc <- function(n,p,conf.level=.95){
alpha <- (1-conf.level)/2
z <- qnorm(1-alpha)
z2 <- z*z
c <- w <- 0
L <- U <- L2 <- U2 <- rep(0,(n+1))

p_w <- 0:n/n  ## 
pt <- p_calc(p_w) ## calculate and store part of prob. for use in function to integrate over m
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt[i] + a*m)*dnorm(m)}
g <- function(m) {(z*sqrt(pnorm(pt[i] + a*m)*(1-pnorm(pt[i] + a*m))/n))*dnorm(m)}
c <- integrate(f,-Inf,Inf)$value   ## center of interval
w <- integrate(g,-Inf,Inf)$value   ## width of interval
L[i] <- c - w
U[i] <- c + w
L2[i] <- c - z*sqrt(c*(1-c)/n)
U2[i] <- c + z*sqrt(c*(1-c)/n)
}
L[which(L<0)] <- 0
L2[which(L2<0)] <- 0
U[which(U>1)] <- 1
U2[which(U2>1)] <- 1

return(c(L,U,L2,U2))
#return(c(cp(n,p,L,U),el(n,p,L,U),cp(n,p,L2,U2),el(n,p,L2,U2)))
}


ac_gofc <- function(n,p,conf.level=.95){
alpha <- (1-conf.level)/2
z <- qnorm(1-alpha)
z2 <- z*z
c <- w <- 0
L <- U <- L2 <- U2 <- rep(0,(n+1))

p_ac <- (((0:n) + 0.5*z2)/(n+z2))  ## add AC adjustments to each possible obs. and form prob.
pt <- p_calc(p_ac) ## calculate and store part of prob. for use in function to integrate over m
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt[i] + a*m)*dnorm(m)}
g <- function(m) {(z*sqrt(pnorm(pt[i] + a*m)*(1-pnorm(pt[i] + a*m))/n))*dnorm(m)}
c <- integrate(f,-Inf,Inf)$value   ## center of interval
w <- integrate(g,-Inf,Inf)$value   ## width of interval
L[i] <- c - w
U[i] <- c + w
L2[i] <- c - z*sqrt(c*(1-c)/n)
U2[i] <- c + z*sqrt(c*(1-c)/n)
}
L[which(L<0)] <- 0
L2[which(L2<0)] <- 0
U[which(U>1)] <- 1
U2[which(U2>1)] <- 1

return(c(L,U,L2,U2))
#return(c(cp(n,p,L,U),el(n,p,L,U),cp(n,p,L2,U2),el(n,p,L2,U2)))
}


bp <- function(n,p) {dbinom(0:n,n,p)}  ## calculate binomial probs.

### cp is function to return coverage prob.
cp <- function(n,p,l,u,probs)	{
pv <- rep(p,length(l))  ## vector of probability
ind <- as.numeric((l<=pv) & (pv<=u))  ## check if p is in interval
return(sum(ind*probs))
#return(sum(ind*gofc_prob_vec(n,p)))
}

### el is function to return expected length
el <- function(n,p,l,u,probs) {sum((u-l)*probs)}

## function calculates part of estimated default prob without M.
p_calc <- function(x) {qnorm(x)*sqrt(1-a^2)} 


### score calculates the score interval

score_gofc <- function(n,p,conf.level=.95){
alpha <- (1-conf.level)/2
z <- qnorm(1-alpha)
z2 <- z*z

c <- w <- 0
L <- U <- L2 <- U2 <- rep(0,(n+1))

pvec <- (0:n)/n
p1 <- (pvec + 0.5 * z2/n)/(1 + z2/n)
pt <- p_calc(p1)
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt[i] + a*m)*dnorm(m)}
g <- function(m) {z * sqrt(((pnorm(pt[i] + a*m) * (1 - pnorm(pt[i] + a*m)) + 0.25 * z2/n)/n)/(1 + z2/n))*dnorm(m)}
c <- integrate(f,-Inf,Inf)$value   ## center of interval
w <- integrate(g,-Inf,Inf)$value
L[i] <- c - w
U[i] <- c + w
L2[i] <- c - (z * sqrt((c * (1 - c) + 0.25 * z2/n)/n)/(1 + z2/n))
U2[i] <- c + (z * sqrt((c * (1 - c) + 0.25 * z2/n)/n)/(1 + z2/n))
	}
L[which(L<0)] <- 0
L2[which(L2<0)] <- 0
U[which(U>1)] <- 1
U2[which(U2>1)] <- 1
return(c(L,U,L2,U2))
#return(c(cp(n,p,L,U),el(n,p,L,U),cp(n,p,L2,U2),el(n,p,L2,U2)))
}


###  mod_score_gofc is version of modified score

mod_score_gofc <- function(n,p,conf.level=.95){
alpha <- (1-conf.level)/2
z <- qnorm(1-alpha)
z2 <- z*z

c <- w <- 0
L <- U <- L2 <- U2 <- rep(0,(n+1))

x <- 0:n
p1 <- rep(0,n+1)
for (j in (x+1)) {
temp <- binomCI(x[j],n,conf.level=conf.level,"modified wilson")$CI
p1[j] <- 0.5*(temp[2]+temp[1])
	}
pt <- p_calc(p1)
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt[i] + a*m)*dnorm(m)}
g <- function(m) {z * sqrt(((pnorm(pt[i] + a*m) * (1 - pnorm(pt[i] + a*m)) + 0.25 * z2/n)/n)/(1 + z2/n))*dnorm(m)}
c <- integrate(f,-Inf,Inf)$value   ## center of interval
w <- integrate(g,-Inf,Inf)$value
L[i] <- c - w
U[i] <- c + w
L2[i] <- c - (z * sqrt((c * (1 - c) + 0.25 * z2/n)/n)/(1 + z2/n))
U2[i] <- c + (z * sqrt((c * (1 - c) + 0.25 * z2/n)/n)/(1 + z2/n))
	}

L[which(L<0)] <- 0
L2[which(L2<0)] <- 0
U[which(U>1)] <- 1
U2[which(U2>1)] <- 1

return(c(L,U,L2,U2))
#return(c(cp(n,p,L,U),el(n,p,L,U),cp(n,p,L2,U2),el(n,p,L2,U2)))
}


### jeff calculates the jefferys interval

jeff_gofc <- function(n,p,conf.level=.95){
L <- U <- rep(0,n+1)
x1 <- (0:n) + 0.5
x2 <- n +0.5 - (0:n)   
pt1 <- qnorm(qbeta((1-conf.level)/2, x1, x2))*sqrt(1-a^2)
pt2 <- qnorm(qbeta(1-(1-conf.level)/2, x1, x2))*sqrt(1-a^2)
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt1[i] + a*m)*dnorm(m)}
g <- function(m) {pnorm(pt2[i] + a*m)*dnorm(m)}
L[i] <- integrate(f,-Inf,Inf)$value   
U[i] <- integrate(g,-Inf,Inf)$value  
}
L[1]<-0
U[n+1]<-1
return(c(L,U))
#return(c(cp(n,p,L,U),el(n,p,L,U)))
}

### jeff calculates the jefferys interval

mod_jeff_gofc <- function(n,p,conf.level=.95){
L <- U <- rep(0,n+1)
x1 <- (0:n) + 0.5
x2 <- n +0.5 - (0:n)   
pt1 <- qnorm(qbeta((1-conf.level)/2, x1, x2))*sqrt(1-a^2)
pt2 <- qnorm(qbeta(1-(1-conf.level)/2, x1, x2))*sqrt(1-a^2)
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt1[i] + a*m)*dnorm(m)}
g <- function(m) {pnorm(pt2[i] + a*m)*dnorm(m)}
L[i] <- integrate(f,-Inf,Inf)$value   
U[i] <- integrate(g,-Inf,Inf)$value  
}
L[1]<-0
U[n+1]<-1
L[2]<-0
U[n]<-1
L[n+1]<-(conf.level/2)^(1/100)
U[1]<-1-(conf.level/2)^(1/100)
return(c(L,U))
#return(c(cp(n,p,L,U),el(n,p,L,U)))
}


### mod_jeff_gofc calculates the modified jefferys interval

mod_jeff_gofc <- function(n,p,conf.level=.95){
L <- U <- rep(0,n+1)
x <- 0:n
pt1 <- pt2 <- rep(0,n+1)
for (j in (x+1)) {
temp <- binomCI(x[j],n,conf.level=conf.level,"modified jeffreys")$CI  
pt1[j] <- qnorm(temp[1])*sqrt(1-a^2)
pt2[j] <- qnorm(temp[2])*sqrt(1-a^2)
}
for (i in 1:(n+1)) {
f <- function(m) {pnorm(pt1[i] + a*m)*dnorm(m)}
g <- function(m) {pnorm(pt2[i] + a*m)*dnorm(m)}
L[i] <- integrate(f,-Inf,Inf)$value   
U[i] <- integrate(g,-Inf,Inf)$value  
}
#L[1]<-0
#U[n+1]<-1
return(c(L,U))
#return(c(cp(n,p,L,U),el(n,p,L,U)))
}

    

### Function zl returns the lower and upper bounds of the Zhou-Li interval

zl_gofc <- function(n,p,conf.level=0.95){
L <- U <- L2 <- U2 <- rep(0,n+1)
pvec <- (0:n)/n
pvec[1] <- 0.5/(n+1)
pvec[n+1] <- (n+0.5)/(n+1)
pt <- p_calc(pvec)
z <- qnorm(1-(1-conf.level)/2)
a<- -1/6
b<- 1/6

for (i in 1:(n+1)) {
ptnorm <- function(m) pnorm(pt[i] + a*m)
gam <- function(m) ifelse(ptnorm(m)==0.5,2,(1-2*ptnorm(m))/sqrt(ptnorm(m)*(1-ptnorm(m))))

ginv <- function(m) return(n^{1/2}*(a*gam(m))^{-1}*(sign(1+3*a*gam(m)*(n^{-1/2}*z-n^{-1}*b*gam(m)))*abs((1+3*a*gam(m)*(n^{-1/2}*z-n^{-1}*b*gam(m))))^{1/3}-1))
ginv2 <- function(m) return(n^{1/2}*(a*gam(m))^{-1}*(sign(1+3*a*gam(m)*(n^{-1/2}*(-z)-n^{-1}*b*gam(m)))*abs((1+3*a*gam(m)*(n^{-1/2}*(-z)-n^{-1}*b*gam(m))))^{1/3}-1))
 
lf <- function(m) dnorm(m)*exp(log(ptnorm(m)/(1-ptnorm(m)))-(n*ptnorm(m)*(1-ptnorm(m)))^{-1/2}*ginv(m))/(1+exp(log(ptnorm(m)/(1-ptnorm(m)))-(n*ptnorm(m)*(1-ptnorm(m)))^{-1/2}*ginv(m)))
L[i] <- integrate(lf,-20,20)$value

uf <- function(m) dnorm(m)*exp(log(ptnorm(m)/(1-ptnorm(m)))-(n*ptnorm(m)*(1-ptnorm(m)))^{-1/2}*ginv2(m))/(1+exp(log(ptnorm(m)/(1-ptnorm(m)))-(n*ptnorm(m)*(1-ptnorm(m)))^{-1/2}*ginv2(m)))
U[i] <- integrate(uf,-20,20)$value

f <- function(m) {pnorm(pt[i] + a*m)*dnorm(m)}
theta <- integrate(f,-Inf,Inf)$value

if (theta==0.5) {gamma<-2} else {gamma<-(1-2*theta)/sqrt(theta*(1-theta))}   ### The value of gamma defined in paper by Zhou-Li is undefined for phat=0.5 but converges to 2 as phat-->0.5, so define to be 2.

ginv3 <- function(x) return(n^{1/2}*(a*gamma)^{-1}*(sign(1+3*a*gamma*(n^{-1/2}*x-n^{-1}*b*gamma))*abs((1+3*a*gamma*(n^{-1/2}*x-n^{-1}*b*gamma)))^{1/3}-1))

L2[i] <- exp(log(theta/(1-theta))-(n*theta*(1-theta))^{-1/2}*ginv3(z))/(1+exp(log(theta/(1-theta))-(n*theta*(1-theta))^{-1/2}*ginv3(z)))
U2[i] <- exp(log(theta/(1-theta))-(n*theta*(1-theta))^{-1/2}*ginv3(-z))/(1+exp(log(theta/(1-theta))-(n*theta*(1-theta))^{-1/2}*ginv3(-z)))

}
return(c(L,U,L2,U2))
#return(c(cp(n,p,L,U),el(n,p,L,U),cp(n,p,L2,U2),el(n,p,L2,U2)))
}

exact_gofc <- function(n,p,conf.level=0.95) {
alpha2 <- 1-(1-conf.level)/2
tol <- 0.00001
L <- U <- rep(0,n+1)
L[1] <- -Inf
U[n+1] <- Inf

cum_prob <- function(x,d){
ptemp <- rep(0,d)
for (j in 0:(d-1)) {
f <- function(m) {pnorm((x-a*m)/sqrt(1-a^2))^j*(1-pnorm((x-a*m)/sqrt(1-a^2)))^(n-j)*dnorm(m)}
ptemp[j+1] <- choose(n,j)* integrate(f,-Inf, Inf)$value
}
return(sum(ptemp))
}

for (i in 1:n) {  
d <- i
x <- -6
y <- 6
theta <- (x+y)/2
iter <- 1
while((abs(cum_prob(theta,d)- alpha2) > tol) && (iter<50)) {
if (cum_prob(theta,d) < alpha2) y <- theta
else x <- theta
theta <- (x+y)/2
iter <- iter +1
}
L[i+1] <- theta

}

cum_prob2 <- function(x,d){
ptemp <- rep(0,d+1)
for (j in 0:d) {
f <- function(m) {pnorm((x-a*m)/sqrt(1-a^2))^j*(1-pnorm((x-a*m)/sqrt(1-a^2)))^(n-j)*dnorm(m)}
ptemp[j+1] <- choose(n,j) * integrate(f,-Inf, Inf)$value
}
return(sum(ptemp))
}

for (i in 0:(n-1)) {  
d2 <- i
x2 <- -6
y2 <- 6
theta2 <- (x2+y2)/2
iter2 <- 1
while((abs(cum_prob2(theta2,d2) - (1-alpha2)) > tol) && (iter2 < 50)) {
if (cum_prob2(theta2,d2) < (1-alpha2)) y2 <- theta2
else x2 <- theta2
theta2 <- (x2+y2)/2
iter2 <- iter2 +1
}
U[i+1] <- theta2

}
L <- pnorm(L)
U <- pnorm(U)
return(c(L,U))
#return(c(cp(n,p,L,U),el(n,p,L,U)))
}