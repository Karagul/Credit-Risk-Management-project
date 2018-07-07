
library(mvtnorm)

############################## inverse survival function; takes #probability as input returns random variate

inv_surv <- function(p) return((-1/lambda)*log(1-p))

copula_sim <- function(n, type, corr, p, dof, runs=10000) {
# n = number of assets
# type = g for gaussian copula, t for t copula
# corr = correlation parameter
# p = probability of default
# dof = degrees of freedom for t copula
# runs = number of simulation runs, default 10000

## generate correlation matrix with equal correlation between assets

corr_mat <- matrix(rep(corr,n^2),nrow=n)
diag(corr_mat) <- rep(1,n)

## depending on input 'g' or 't', generate multivariate random variates
## from appropriate distribution with above correlation matrix.
## for each dimension, check if this would 'default' by comparing to
## inverse of 'p' value supplied by user, appropriate univariate distribution.
## for each value less than this value, assign 1 for default.
## do this 'runs' time, default being 10000, and store in matrix results
 
if (type == 'g') 
 results <- matrix(as.numeric(rmvnorm(runs,sigma=corr_mat,method="chol") < qnorm(p)),ncol=n)
else if (type == 't')
  results <- matrix(as.numeric(rmvt(runs,sigma=corr_mat,df=dof) < qnorm(p)),ncol=n)
else
return(paste("please specify a either 'g' for gaussian copula or 't' for t"))
 
return(sum(rowSums(results))/runs)
}
  
count<-function(x) {
y<-vector("numeric")
for (i in 1:length(t)) y[i]<- length(x[x<g_def_prob(t[i]]))
return(y)
}
  
	 
		

