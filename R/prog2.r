

## rm(list=ls(all=TRUE))

## input the correlation between bond price movements, mean and standard deviation of recovery rate for defaulted bonds,
## and S&P ratings

corr <- .20
recov_rate <- .5113
recov_sd <- .2545
ratings <- c("AAA","AA","A","BBB","BB","B","CCC","D")

## input the transition matrix

trans_mat <- matrix(c(90.81,0.70,0.09,0.02,0.03,0,0.22,
8.33,90.65,2.27,0.33,0.14,0.11,0,
0.68,7.79,91.05,5.95,0.67,0.24,0.22,
0.06,0.64,5.52,86.93,7.73,0.43,1.30,
0.12,0.06,0.74,5.30,80.53,6.48,2.38,
0,0.14,0.26,1.17,8.84,83.46,11.24,
0,0.02,0.01,0.12,1,4.07,64.86,
0,0,0.06,0.18,1.06,5.2,19.79),
nrow=7,ncol=8,
dimnames=list(ratings[1:7],ratings))

## put trans_mat in percentage format

trans_mat<-trans_mat/100

## input the forward-zero curve

zero_curve <- matrix(c(.036,.0365,.0372,.041,.0555,.0605,.1505,
.0417,.0422,.0432,.0467,.0602,.0702,.1502,
.0473,.0478,.0493,.0525,.0678,.0803,.1403,
.0512,.0517,.0532,.0563,.0727,.0852,.1352),
nrow=7,ncol=4,
dimnames=list(ratings[1:7],
c("yr1","yr2","yr3","yr4")))

## input the bonds in portfolio, the create a portfolio to hold individual bonds

ibm <- list(sym='ibm',face=1000,coupon=.05,rating='A',mat=3,quant=1,stdv=0.015498)
ko <- list(sym='ko',face=1000, coupon=.06,rating='BBB',mat=4,quant=2,stdv=0.013581)
f <- list(sym='f',face=1000,coupon=.07,rating='B',mat=5,quant=3,stdv=0.038252)
portfolio <- list(ibm,ko,f)
names(portfolio) <- c('ibm','ko','f')

## function to caclulate the present value of a bond given a rating.  input is bond and rating, output is present value

pres_val <- function(bond,rating) {
	v <- bond$face*bond$coupon
	for (i in 1:(bond$mat-2)) {
		v <- v + (bond$face*bond$coupon)/((1 + zero_curve[rating,i])^i)
	}
	v <- v + (bond$face + (bond$face*bond$coupon))/((1 + zero_curve[rating,bond$mat-1])^(bond$mat-1))
	return(v)
}

## create a vector of ratings in porfolio, then extract unique ratings

bond_ratings<-vector(length=length(portfolio))
for (i in 1:length(portfolio)) {
bond_ratings[i]<-portfolio[[i]]$rating
}
uniq_rat<-unique(bond_ratings)

## for each unique rating in portfolio, calculate value of rating change from initial rating to each rating in rating system
## e.g. for bond initially in AAA, calc value of bond for change to AAA, AA, A...
## store in matrix value_dsn

value_dsn <- matrix(rep(0,length(uniq_rat)*length(ratings)),nrow=length(uniq_rat),
ncol=length(ratings),dimnames=list(uniq_rat,ratings))


for (i in 1:length(portfolio)) {
	for (j in 1:(length(ratings)-1)) {
		value_dsn[i,j] <- pres_val(portfolio[[i]],ratings[j])
	}
	value_dsn[i,length(ratings)] <- recov_rate*portfolio[[i]]$face
}

## calculate the threshold values for rating changes for each bond.  thresh1 is with standard deviation 1 and 
## thresh is thresh1 multiplied by bond stand. dev.  
## fix: this should be done for unique ratings, not every bond in portfolio

thresh1<-vector("list", length=length(portfolio))
thresh<-vector("list", length=length(portfolio))
for (i in 1:length(portfolio)) {
	temp<-vector(length=length(ratings)+1)
	temp[1]<-Inf
	temp[length(ratings)+1]<--Inf
	for (j in 2:(length(ratings))) {
		temp[j]<-qnorm(1 - sum(trans_mat[portfolio[[i]]$rating,1:(j-1)]))
	}
	names(temp)<-c(ratings,"subD")
	thresh1[[i]]<-temp
}
names(thresh1)<-names(portfolio)

for (i in 1:length(portfolio)) {
thresh[[i]]<-thresh1[[i]]*portfolio[[i]]$stdv
}

## calculate transition matrix for each two bond combination in portfolio.
## prob_mat is a list containing each transition matrix

prob_mat<-vector("list", length=choose(length(portfolio),2))
count<-1
prob_mat_names<-vector(length=choose(length(portfolio),2))

for (i in 1:((length(portfolio)-1))) {
	for (j in (i+1):(length(portfolio))) {
		covar<-matrix(c((portfolio[[i]]$stdv)^2,corr*portfolio[[i]]$stdv*portfolio[[j]]$stdv,
			corr*portfolio[[i]]$stdv*portfolio[[j]]$stdv,(portfolio[[j]]$stdv)^2),nrow=2)
		prob_mat_names[count]<-paste(names(portfolio[i]),names(portfolio[j]),sep="_")
		temp<-matrix(rep(0,64), nrow=8)
		for (k in 1:(length(ratings))) {
			for (l in 1:(length(ratings))) {
				temp[k,l]<-pmvnorm(lower=c(thresh[[i]][k+1],thresh[[j]][l+1]),
						upper=c(thresh[[i]][k],thresh[[j]][l]),sigma=covar)
			}
		}
		prob_mat[[count]]<-temp
		dimnames(prob_mat[[count]])<-list(ratings,ratings)
		count<-count+1
	}
}

names(prob_mat)<-prob_mat_names

## function mean1 takes as input a bond and returns expected value

mean1 <- function(bond) {
	mean<-0
	for (i in 1:length(ratings)) {
	mean <- mean + (value_dsn[bond$rating,i] * trans_mat[bond$rating,i])
	}
return(mean)
}

## function mean2 takes as input 2 bonds and returns expected value of the 2 bond portfolio

mean2 <- function(bond1,bond2) {
	mean<-0
	for (i in 1:length(ratings)) {
		for (j in 1:length(ratings)) {
			mean <- mean + (value_dsn[bond1$rating,i] + value_dsn[bond2$rating,j]) * 
				prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][i,j]
		}
	}
return(mean)
}

## function mean1_sq takes as input a bond and returns expected squared value
## ## note here the standard deviation of recovery rate is included in default scenarios to capture variability of recovery

mean1_sq <- function(bond) {
	mean<-0
	for (i in 1:(length(ratings)-1)) {
	mean <- mean + (((value_dsn[bond$rating,i])^2) * trans_mat[bond$rating,i])
	}
	mean <- mean + (value_dsn[bond$rating,length(ratings)]^2 + (bond$face*recov_sd)^2) * trans_mat[bond$rating,length(ratings)]
return(mean)
}

## function mean2_sq takes as input 2 bonds and returns expected squared value of 2 bond portfolio
## note here the standard deviation of recovery rate is included in default scenarios to capture variability of recovery

mean2_sq <- function(bond1,bond2) {
        mean<-0
        for (i in 1:(length(ratings)-1)) {
                for (j in 1:(length(ratings)-1)) {
                        mean <- mean + (value_dsn[bond1$rating,i] + value_dsn[bond2$rating,j])^2 * 
                                prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][i,j]
                }
        }
	for (k in 1:(length(ratings)-1)) {
		mean <- mean + ((value_dsn[bond1$rating,length(ratings)] + value_dsn[bond2$rating,k])^2 + (recov_sd*bond1$face)^2) *
			prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][length(ratings),k]
	}
	for (l in 1:(length(ratings)-1)) {
		mean <- mean + ((value_dsn[bond1$rating,l] + value_dsn[bond2$rating,length(ratings)])^2 + (recov_sd*bond2$face)^2) *
			prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][l,length(ratings)]
	}
	mean <- mean + ((recov_sd*bond1$face)^2 + (recov_sd*bond2$face)^2 + (value_dsn[bond1$rating,length(ratings)]
	 + value_dsn[bond2$rating,length(ratings)])^2) * prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][length(ratings),length(ratings)]
	
return(mean)
}

## function sd1 takes input bond and returns standard deviation of value of the bond

sd1 <- function(bond) {
	sd <- sqrt(mean1_sq(bond) - (mean1(bond))^2)	
	return(sd)
	}

## function sd2 takes as input 2 bonds and returns standard deviation of value of 2 bond portfolio

sd2 <- function(bond1,bond2) {
	sd <- sqrt(mean2_sq(bond1,bond2) - (mean2(bond1,bond2))^2)	
	return(sd)
	}

## function port_sd takes as input a portfolio of bonds and returns standard deviation of portfolio value

port_sd <- function(portfolio) {
	dev1<-0
	dev2<-0
	for (i in 1:((length(portfolio)-1))) {
		for (j in (i+1):(length(portfolio))) {
			dev2 <- dev2 + sd2(portfolio[[i]],portfolio[[j]])
		}
	}
	for (k in 1:length(portfolio)) {
		dev1 <- dev1 + sd1(portfolio[[k]])
	}
return(dev2-dev1)
}

## analytic_sd stores value of portfolio standard deviation

analytic_sd <- port_sd(portfolio)



## BEGIN SIMULATION CODE


## input correlation matrix
	
corr_mat <- matrix(rep(corr,(length(portfolio)^2)),nrow=length(portfolio))
diag(corr_mat) <- rep(1,length(portfolio))

## function sim_rat takes as input a value (here a standard normal random variate representing an asset return)  
## and an integer representing a bond in the portfolio and returns what rating results from the value

sim_rat <- function(value,z) {
	k=2
	while (k <= (length(ratings)+1)) {
		if (value > thresh1[[z]][k]) {
			rate=names(thresh1[[z]][k-1])
			break
		}
		k<-k+1
	}
return(rate)
}

## define the number of simulations to run and generate multivariate normal random numbers		

n=100000
sim <- rmvnorm(n, sigma=corr_mat)

## create matrix of simulated portfolio values
## note: if a bond defaults, a random beta variate is generated and used as the recovery rate
## the beta parameters were selected by mean and stand. dev. listed in technical document

sim_val <- matrix(rep(0,n*length(portfolio)), nrow=n, dimnames=list(NULL,names(portfolio)))

for (i in 1:n) {
	for (j in 1:length(portfolio)) {
		if (sim_rat(sim[i,j],j) == 'D') sim_val[i,j] <- rbeta(1,.68,.603)*portfolio[[j]]$face
		else sim_val[i,j] <- value_dsn[portfolio[[j]]$rating,sim_rat(sim[i,j],j)]
	}
}

## sum bond values to get portfolio value for each scenario and add to value matrix

sim_port_val <- vector(length=n)
for (i in 1:n) {
	sim_port_val[i] <- sum(sim_val[i,])
}

sim_val <- cbind(sim_val,sim_port_val)

## sim_quantile function takes as input a percent value and returns the corresponding quantile of simulated distribution

sim_quantile <- function(q) {
	return(quantile(sim_port_val,q,type=8))
}

## sim_sd calculates the simulated portfolio distribution standard deviation

sim_sd <- sd(sim_port_val)

