
mc2 <- function(al,ah,runs=10000) {
corr_mat <- matrix(rep(corr,(N_c)^2),nrow=N_c)
diag(corr_mat) <- rep(1,N_c)
s <- vector("numeric", length=runs)  ## store spread from each run
x<-rep(0,length(t))
nom <- vector("numeric",length=length(t)) ## store nominal value for each time t_i for each run
for (i in 1:runs) {
d <- count3(dt(pnorm(rmvnorm(1,sigma=corr_mat,method="chol"))))
x <- x + d 
#s[i] <- spread(nom)
 }
for (j in 1:(N_c+1)) nv[j] <- rn(j-1,al,ah)
x <- x/runs
for (k in 1:length(d)) nom[k] <- rn(x[k],al,ah)
#return(spread(nom))
return(x)
}


count3 <- function(x) {
y <- vector("numeric")
for (i in 1:length(t)) y[i] <- length(x[x < t[i]])
return(y)
}

dt <- function(p) return(-(1/lambda)*log(1-p))