############  Implied Copula Implementation #############

lam <- vector("numeric")  ## lambda values
lam[1] <- 0.0001
for (index in 1:99) lam[index+1] <- exp(log(lam[index]) + log(1.125) +(log(0.5)-log(1.125))/99)
tr <- c(0,.03,.07,.1,.15,.3)  # vector used for tranche values
tr_sp <- c(.05,.0127,.00355,.00205,.00095,.005) # market quotes for first 5 tranches plus index
tr1_upfront <- .4  # upfront payment on first tranche as percent of principal
start <- rep(0.01,length(lam))  ## initial value
V <- Vcalc()  ## calculate V as a global variable

### Function pp takes k, lambda, and t and returns prob of k defaults by time t for poisson process with rate lambda
pp <- function(k,lam,t) return((exp(-lam*t)*(lam*t)^k)/factorial(k))

### Function nv takes tranche attach points and returns vector of remaining notional after each possible number of defaults
nv <- function(al,ah) {
n <- vector("numeric")
for (i in 1:(N_c+1)) n[i] <- rn(i-1,al,ah)
return(n)
}

### Function val takes vector of nominal values for each time t and market spread and returns contract value
val <- function(nominal,q) {
nominal2 <- nominal[-1]
a <- sum(diff(t)*nominal2*exp(-r*t2))
b <- 0.5*(sum(diff(t)*(-diff(nominal))*exp(-0.5*r*t3)))
c <- sum((-diff(nominal))*exp(-0.5*r*t3))
s <- q*a +q*b - c
return(s)
}

### Function en takes hazard rate lam and tranche points and returns vector of Expected Notional at each time t_i
en <- function(lam,al,ah) {
a <- vector("numeric")
n <- nv(al,ah)
for (i in 1:length(t)) a[i] <- sum(pp(0:N_c,lam,t[i])*n)
return(a)
}

### Function Vcalc calculates matrix of values for each tranche and hazard rate
Vcalc <- function() {
V <- matrix(rep(0,length(lam)*(length(tr))),ncol=length(lam))
for (i in 1:(length(tr)-1)) {
 for (j in 1:length(lam)) V[i,j] <- val(en(lam[j],tr[i],tr[i+1]),tr_sp[i])
 }
for (k in 1:length(lam)) V[length(tr),k] <- val(en(lam[k],0,1),tr_sp[length(tr_sp)])  ## calc the index tranche last
V[1,] <- V[1,] + N_c*K*(tr[2]-tr[1])*tr1_upfront  ## add the upfront payment to 1st tranche
return(V)
}



#########  The count function below is used to take a simulation of default times and count the number of defaults by time t_i
count2 <- function(x) {
y <- vector("numeric")
for (i in 1:length(t)) y[i] <- length(x[x < t[i]])
return(y)
}



## function Vsim simulates contract value for each tranche and hazard by simulating poisson process with correct hazard rate
Vsim <- function(runs=10000) {
V <- matrix(rep(0,length(lam)*length(tr)),ncol=length(lam))
for (i in 1:(length(tr)-1)) {
 for (j in 1:length(lam)) {
  W <- Psim(lam[j],runs)
  V[i,j] <- val(W%*%nv(tr[i],tr[i+1]),tr_sp[i])
  }
 }
for (k in 1:length(lam)) V[length(tr),k] <- val(Psim(lam[k],runs)%*%nv(0,1),tr_sp[length(tr_sp)])
V[1,] <- V[1,] + N_c*K*(tr[2]-tr[1])*tr1_upfront  ## add the upfront payment to 1st tranche 
return(V)
}

### Function Psim takes hazard rate lam and runs and returns matrix of default probs by each time t_i
Psim <- function(lam,runs) {
W <- matrix(rep(0,length(t)*(N_c+1)),nrow=length(t))
  for (i in 1:runs) {
  x <- rexp(N_c,lam)  #simulate N_c 'interarrival' times of default
  y <- cumsum(x) #create vector of actual default times
  d <- count2(y) #create vector of number of defaults by each time t_i
  for (j in 1:length(t)) W[j,d[j]+1] <- W[j,d[j]+1] + 1
  }
W <- W/runs
return(W)
}

fm <- function(x) {
y<-vector("numeric")
for (i in 1:(length(x)-2)) y[i] <- ((x[i+2] - 2*x[i+1] + x[i])^2)/(lam[i+2] - lam[i])
z <- (V%*%x)^2
0.1*sum(y) + sum(z)
}


###Function ec defines the equality constraint that the parameters must sum to 1
ec <- function(x) {
z1 <- sum(x) - 1
return(z1)
}

##Function mopt inputs an initial guess and proceeds with minimization of objective function fm subject to constraints
mopt <- function(init=start) {
UB <- rep(1,length(lam))  ##upper bound on parameters
LB <- rep(0, length(lam)) ##lower bound on parameters
x <- solnp(start,fun=fm,eqfun=ec,eqB=0,LB=LB,UB=UB)
return(x$par)
}





