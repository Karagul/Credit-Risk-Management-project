# local({pkg <- select.list(sort(.packages(all.available = TRUE)),graphics=TRUE)
# + if(nchar(pkg)) library(pkg, character.only=TRUE)})

rm(list=ls(all=TRUE))

corr <- .20
recov_rate <- .5113
recov_sd <- .2545
ratings <- c("AAA","AA","A","BBB","BB","B","CCC","D")

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

trans_mat<-trans_mat/100

zero_curve <- matrix(c(.036,.0365,.0372,.041,.0555,.0605,.1505,
.0417,.0422,.0432,.0467,.0602,.0702,.1502,
.0473,.0478,.0493,.0525,.0678,.0803,.1403,
.0512,.0517,.0532,.0563,.0727,.0852,.1352),
nrow=7,ncol=4,
dimnames=list(ratings[1:7],
c("yr1","yr2","yr3","yr4")))

ibm <- list(sym='ibm',face=1000,coupon=.05,rating='A',mat=3,quant=1,stdv=0.015498)
ko <- list(sym='ko',face=1000, coupon=.06,rating='BBB',mat=4,quant=2,stdv=0.013581)
f <- list(sym='f',face=1000,coupon=.07,rating='B',mat=5,quant=3,stdv=0.038252)
ge <- list(sym='ge',face=100000,coupon=.02,rating='CCC',mat=3,stdv=1)
ge2 <- list(sym='ge2',face=1000,coupon=.02,rating='CCC',mat=3,stdv=0.02)
ge3 <- list(sym='ge3',face=1000,coupon=.02,rating='CCC',mat=4,stdv=0.02)
portfolio <- list(ibm,ko,f,ge,ge2)
names(portfolio) <- c('ibm','ko','f','ge','ge2')

pres_val <- function(bond,rating) {
	if (bond$mat == 1) {
	v <- bond$face*bond$coupon + bond$face
	}
	else if (bond$mat == 2) {
	v <- bond$face*bond$coupon + (bond$face + bond$face*bond$coupon)/zero_curve[rating,1]
	}
	else {
	v <- bond$face*bond$coupon
	for (i in 1:(bond$mat-2)) {
		v <- v + (bond$face*bond$coupon)/((1 + zero_curve[rating,i])^i)
	}
	v <- v + (bond$face + (bond$face*bond$coupon))/((1 + zero_curve[rating,bond$mat-1])^(bond$mat-1))
	}
	return(v)
}

bond_ratings<-vector(length=length(portfolio))
for (i in 1:length(portfolio)) {
bond_ratings[i]<-portfolio[[i]]$rating
}


value_dsn <- matrix(rep(0,length(portfolio)*length(ratings)),nrow=length(portfolio),
ncol=length(ratings),dimnames=list(names(portfolio),ratings))


for (i in 1:length(portfolio)) {
	for (j in 1:(length(ratings)-1)) {
		value_dsn[i,j] <- pres_val(portfolio[[i]],ratings[j])
	}
	value_dsn[i,length(ratings)] <- recov_rate*portfolio[[i]]$face
}


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

mean1 <- function(bond) {
	mean<-0
	for (i in 1:length(ratings)) {
	mean <- mean + (value_dsn[bond$sym,i] * trans_mat[bond$rating,i])
	}
return(mean)
}

mean2 <- function(bond1,bond2) {
	mean<-0
	for (i in 1:length(ratings)) {
		for (j in 1:length(ratings)) {
			mean <- mean + (value_dsn[bond1$sym,i] + value_dsn[bond2$sym,j]) * 
				prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][i,j]
		}
	}
return(mean)
}

mean1_sq <- function(bond) {
	mean<-0
	for (i in 1:(length(ratings)-1)) {
	mean <- mean + (((value_dsn[bond$sym,i])^2) * trans_mat[bond$rating,i])
	}
	mean <- mean + (value_dsn[bond$sym,length(ratings)]^2 + (bond$face*recov_sd)^2) * trans_mat[bond$rating,length(ratings)]
return(mean)
}


mean2_sq <- function(bond1,bond2) {
        mean<-0
        for (i in 1:(length(ratings)-1)) {
                for (j in 1:(length(ratings)-1)) {
                        mean <- mean + (value_dsn[bond1$sym,i] + value_dsn[bond2$sym,j])^2 * 
                                prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][i,j]
                }
        }
	for (k in 1:(length(ratings)-1)) {
		mean <- mean + ((value_dsn[bond1$sym,length(ratings)] + value_dsn[bond2$sym,k])^2 + (recov_sd*bond1$face)^2) *
			prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][length(ratings),k]
	}
	for (l in 1:(length(ratings)-1)) {
		mean <- mean + ((value_dsn[bond1$sym,l] + value_dsn[bond2$sym,length(ratings)])^2 + (recov_sd*bond2$face)^2) *
			prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][l,length(ratings)]
	}
	mean <- mean + ((recov_sd*bond1$face)^2 + (recov_sd*bond2$face)^2 + (value_dsn[bond1$sym,length(ratings)]
	 + value_dsn[bond2$sym,length(ratings)])^2) * prob_mat[paste(bond1$sym,bond2$sym,sep='_')][[1]][length(ratings),length(ratings)]
	
return(mean)
}

sd1 <- function(bond) {
	sd <- sqrt(mean1_sq(bond) - (mean1(bond))^2)	
	return(sd)
	}

sd2 <- function(bond1,bond2) {
	sd <- sqrt(mean2_sq(bond1,bond2) - (mean2(bond1,bond2))^2)	
	return(sd)
	}

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

analytic_sd <- port_sd(portfolio)  
	
corr_mat <- matrix(rep(corr,(length(portfolio)^2)),nrow=length(portfolio))
diag(corr_mat) <- rep(1,length(portfolio))

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
		

n=20000
sim <- rmvnorm(n, sigma=corr_mat)

sim_val <- matrix(rep(0,n*length(portfolio)), nrow=n, dimnames=list(NULL,names(portfolio)))

for (i in 1:n) {
	for (j in 1:length(portfolio)) {
		if (sim_rat(sim[i,j],j) == 'D') sim_val[i,j] <- rbeta(1,.68,.603)*portfolio[[j]]$face
		else sim_val[i,j] <- value_dsn[portfolio[[j]]$sym,sim_rat(sim[i,j],j)]
	}
}
sim_port_val <- vector(length=n)
for (i in 1:n) {
	sim_port_val[i] <- sum(sim_val[i,])
}



sim_val <- cbind(sim_val,sim_port_val)

sim_quantile <- function(q) {
	return(quantile(sim_port_val,q,type=8))
}
sim_sd <- sd(sim_port_val)

