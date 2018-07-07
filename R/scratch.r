nom_val2 <- function(al,ah) {
d <- vector("numeric")
nv <- vector("numeric") 
for (i in 1:(N_c+1)) nv[i] <- rn(i-1,al,ah)
for (j in 1:length(t)) nom[j] <- sum(nv*rFLD(t[j]))
return(nom)
}

mc4 <- function(al,ah,runs=10000) {
corr_mat <- matrix(rep(corr,(N_c)^2),nrow=N_c)
diag(corr_mat) <- rep(1,N_c)
W <- matrix(rep(0,length(t)*(N_c+1)),nrow=length(t))
x <- rep(0,length(t))
s <- vector("numeric", length=runs)  ## store spread from each run
nom <- vector("numeric",length=length(t)) ## store nominal value for each time t_i for each run
for (i in 1:runs) {
d <- count(rmvnorm(1,sigma=corr_mat,method="chol"))
for (j in 1:length(t)) W[j,d[j]+1] <- W[j,d[j]+1] + 1
 }
W <- W/runs
for (i in 1:(N_c+1)) nv[i] <- rn(i-1,al,ah)
nom <- W%*%nv
return(nom)
}