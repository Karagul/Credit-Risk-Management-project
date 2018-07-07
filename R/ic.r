############  Implied Copula Implementation #############

lam <- vector("numeric")  ## lambda values
lam[1] <- 0.001
for (index in 1:99) lam[index+1] <- exp(log(lam[index]) + log(1.125) +(log(0.5)-log(1.125))/99)
tr <- c(0,.03,.07,.1,.15,.3,1.3)  # vector used for tranche values
tr_sp <- c(.05,.0127,.00355,.00205,.00095,.005) # market quotes for first 5 tranches plus index
C <- diag(length(lam))  ## matrix of constraints on probabilities
C <- rbind(C,matrix(c(rep(1,length(lam)),rep(-1,length(lam))),nrow=2,byrow=TRUE))
con <- c(rep(0,length(lam)),1,-1.0011)  ## vector used for constraints
init <- rep(0.01,length(lam))  ## initial value

### Function pp takes k, lambda, and t and returns prob of k defaults by time t for poisson process with rate lambda
pp <- function(k,lam,t) return((exp(-lam*t)*(lam*t)^k)/factorial(k))

### Function nv takes tranche attach points and returns vector of remaining notional after each possible number of defaults
nv <- function(al,ah) {
n <- vector("numeric")
for (i in 1:(N_c+1)) n[i] <- rn(i-1,al,ah)
return(n)
}

### Function nominal takes hazard rate lam and tranche attach points and returns expected nominal value for each time t_i
nominal <- function(lam,al,ah) {
nom <- vector("numeric")
for (i in 1:length(t)) nom[i] <- sum(pp(0:125,lam,t[i])*nv(al,ah))
return(nom)
}

### Function val takes vector of nominal values for each time t and market spread and returns contract value
val <- function(nominal,q) {
nominal2 <- nominal[-1]
a <- sum(diff(t)*nominal2*exp(-r*t2))
b <- 0.5*(sum(diff(t)*(-diff(nominal))*exp(-0.5*r*t3)))
c <- (1-R)*sum((-diff(nominal))*exp(-0.5*r*t3))
s <- q*a +q*b - c
return(s)
}



#########  The count function below is used to take a simulation of default times and count the number of defaults by time t_i
count2 <- function(x) {
y <- vector("numeric")
for (i in 1:length(t)) y[i] <- length(x[x < t[i]])
return(y)
}

## function spread takes nominal values at each time and spread quote and calculates PV of each leg and returns contract value
spread2 <- function(nominal,q) {
nominal2 <- nominal[-1]
a <- sum(diff(t)*nominal2*exp(-r*t2))
b <- 0.5*(sum(diff(t)*(-diff(nominal))*exp(-0.5*r*t3)))
c <- (1-R)*sum((-diff(nominal))*exp(-0.5*r*t3))
s <- q*a +q*b - c
return(s)
}


## function v takes number of runs,hazard rate,tranche attach points and tranche spread and simulates contract value
v <- function(runs,haz,ap,dp,sp) {
s <- vector("numeric")
for (k in 1:runs) {
 nom <- vector("numeric")
 x <- rexp(125,haz)  #simulate 125 'interarrival' times of default
 y <- cumsum(x) #create vector of actual default times
 z <- count2(y) #create vector of number of defaults by each time t_i
 for (l in 1:length(z)) nom[l] <- rn(z[l],ap,dp)
 s[k] <- spread2(nom,sp)
 }
spread <- sum(s)/runs
return(spread)
}


# First build a matrix of values for each tranche and each lambda value
#V <- matrix(rep(0,length(lam)*(length(tr)-1)),ncol=length(lam))
#for (i in 1:(length(tr)-1)) {
# for (j in 1:length(lam)) V[i,j] <- v(100,lam[j],tr[i],tr[i+1],tr_sp[i])
#}

fm <- function(x) {
y<-vector("numeric")
for (i in 1:(length(x)-2)) y[i] <- ((x[i+2] - 2*x[i+1] + x[i])^2)/(lam[i+2] - lam[i])
z <- (V%*%x)^2
.01*sum(y) + sum(z)
}

grad<-function(x) {
g<-vector("numeric")
g[1] <- (1/(lam[3]-lam[1]))*(2*(x[3]-2*x[2]+x[1]))
g[2] <- (-4/(lam[3]-lam[1]))*(x[3]-2*x[2]+x[1]) + (2*(lam[4]-lam[2]))*(x[4]-2*x[3]+x[2])
for (i in 3:(length(x)-2)) g[i] <- (1/(lam[i]-lam[i-2]))*(2*(x[i]-2*x[i-1]+x[i-2])) 
 +(-4/(lam[i+1]-lam[i-1]))*(x[i+1]-2*x[i]+x[i-1]) + (2*(lam[i+2]-lam[i]))*(x[i+2]-2*x[i+1]+x[1])
g[length(x)-1]<- (-4/(lam[length(x)]-lam[length(x)-2]))*(x[length(x)]-2*x[length(x)-1]+x[length(x)-2])
 + (2*(lam[length(x)-1]-lam[length(x)-3]))*(x[length(x)-1]-2*x[length(x)-2]+x[length(x)-3])
g[length(x)] <- (1/(lam[length(x)]-lam[length(x)-2]))*(2*(x[length(x)]-2*x[length(x)-1]+x[length(x)-2]))
g
}

ww<-constrOptim(theta=init,f=fm,grad=grad,ui=C,ci=con,method="BFGS")






