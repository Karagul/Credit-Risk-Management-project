

## change mc to be able to calc diff tranches without re-running sims

N_c <- 50;                ## number of credits in portfolio
K <- 1			   ## notional of each credit
lambda <- .05;             ## parameter for exponential dist. used for default probability
		               ## could be vectorized if each credit has own parameter
T <- 5;	               ## Time until maturity of portfolio in years
t <- seq(0,T,.25)	   ## divide T into quarters for payment calculation
t1 <- t[-length(t)]   ## t1,t2,t3 will be used in spread calc.
t2 <- t[-1]
t3 <- t2 + t1
R <- 0.4	               ## recovery rate
corr <- .2       	         ## correlation between credits
z <- c(seq(-4,-2,.5),seq(-1.9,1.9,.1),seq(2,4,.5))         ## vector of market factor values
z1 <- z[-length(z)]        ## sequence used in integrating over z
z2 <- z[-1]                ## same as above
z3 <- pnorm(z2)-pnorm(z1)  ## same as above
r <- .0342			   ## risk free interest rate                



## given time t, return prob of default by t
def_prob <- function(t) return(1 - exp(-lambda*t)); ## given time t, return prob
								    ## of default by t
## map default prob to gaussian dist.
g_def_prob <- function(t) return(qnorm(def_prob(t))); 

## function dt takes probability and calculates default time by inverting expo dist.
dt <- function(p) return(-(1/lambda)*log(1-p))

## calculate conditional default probability for a credit in portfolio
cond_def_prob <- function(time, factor) return(pnorm((g_def_prob(time) - corr*factor)/(sqrt(1-corr^2))));

########  Calculate conditional expected defaults using recursive method ########
## function rCLD calculates Conditional Loss Distribution (CLD) recursivley
## given time and factor values

rCLD <- function(time,factor) {
p <- cond_def_prob(time,factor)
## initialize values; f[i] will hold prob. of f(i-1) defaults
f <- vector("numeric", length=N_c+1); 
f[1]<-1;
for (i in 2:(N_c+1)) f[i] <- 0;

## now recursively calculate cond loss by adding one credit at a time
for (j in 2:(N_c+1)) {
 for (k in j:2) {f[k] <- f[k]*(1-p) + f[k-1]*p}
 f[1] <- f[1]*(1-p)
 }
return(f)
}

########  Calculate total expected defaults using recursive method by integrating over market factor########
## function rFLD takes as input time T and and returns the expected number of defaults by time T

rFLD <- function(time) {
d <- mapply(rCLD,time=time,factor=z2) ## calc. prob. of number of defaults for each factor in z2
return(d%*%z3) ## calculate and return full dist. of defaults by numerically integrating over market factor
}

########  Calculate remaining notional after j defaults ########
## input number defaults j,low attach point, al, and high attach point ah

rn <- function(j,al,ah) {
nl <- al*N_c/(1-R)
nh <- ah*N_c/(1-R)
if (j < ceiling(nl)) p <- (ah-al)*N_c*K
 else if (ceiling(nl) <= j & j < ceiling(nh)) p <- ah*N_c*K - j*(1-R)*K
 else if (j >= ceiling(nh)) p <- 0
return(p)
}

##  function nom_val takes trance points and returns expected nominal value at each time of t
nom_val <- function(al,ah) {
nom <- vector("numeric", length=length(t))
nv <- vector("numeric", length=N_c+1)
for (j in 1:(N_c+1)) nv[j] <- rn(j-1,al,ah)
for (i in 1:length(t)) nom[i] <- sum(nv*rFLD(t[i]))
return(nom)
}

## function spread takes nominal values at each time and calculates PV of each leg and returns spread
spread <- function(nominal) {
nominal2 <- nominal[-1]
a <- sum(diff(t)*nominal2*exp(-r*t2))
b <- 0.5*(sum(diff(t)*(-diff(nominal))*exp(-0.5*r*t3)))
c <- (1-R)*sum((-diff(nominal))*exp(-0.5*r*t3))
s <- c/(a+b)
return(s)
}

############ function recur takes tranche points and returns spread

recur <- function(al,ah) return(spread(nom_val(al,ah)))


########  Calculate total expected defaults using monte carlo simulation ########

mc2 <- function(al,ah,runs=1000) {
s <- vector("numeric", length=runs)  ## store spread from each run
nom <- vector("numeric",length=length(t)) ## store nominal value for each time t_i for each run
x <- vector("numeric")
d <- vector("numeric")
for (i in 1:runs) {
M <- rnorm(1)
x <- corr*M + (sqrt(1-corr^2))*rnorm(N_c)
d <- count3(dt(pnorm(x)))
for (k in 1:length(d)) nom[k] <- rn(d[k],al,ah) 
s[i] <- spread(nom)
 }
return(sum(s)/runs)
}

#########  The count function below is used to take a monte carlo run and count the number of defaults by time t_i
count <- function(x) {
y <- vector("numeric")
for (i in 1:length(t)) y[i] <- length(x[x < g_def_prob(t[i])])
return(y)
}

#####  count3 function takes vector of default times and counts how many defaults have occured by each t_i
count3 <- function(x) {
y <- vector("numeric")
for (i in 1:length(t)) y[i] <- length(x[x < t[i]])
return(y)
}






