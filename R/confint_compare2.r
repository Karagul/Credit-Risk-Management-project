
### function confint_compare2 inputs probability of defaults for two classes (p1 and p2) and number of observations (n1 and n2) along with a confidence level.  
### Agresti-Coull, Wald, Score, Jefferys, and Zhou-Li two-sided confidence intervals are formed after constrained estimation of p1 and p2.  If p2 > p1, then we set p1=p2
### The output is a matrix of coverage probability and expected length for each method

confint_compare2 <- function(p1,n1,p2,n2,conf.level){

alpha <- 1 - conf.level;
alpha2 <- 0.5 * alpha;

ac_cp1<-ac_cp2<-w_cp1<-w_cp2<-s_cp1<-s_cp2<-j_cp1<-j_cp2<-zl_cp1<-zl_cp2<-0;  ## set all coverage probabilities to zero
ac_el1<-ac_el2<-w_el1<-w_el2<-s_el1<-s_el2<-j_el1<-j_el2<-zl_el1<-zl_el2<-0;  ## set all expected length to zero

for (x2 in 0:n2){											## Loop through each value of the number of defaults that are possible
  for (x1 in 0:n1) {
    if (x1/n1 <= x2/n2) {X1<-x1;X2<-x2;N1<-n1;N2<-n2} 				## Check to see if constraints on esimates are satisfied; if not -->make both equal to higher estimate
	  else {X1<-X2<-x2; N1<-N2<-n2}

### Given the estimates above, the following code computes the various intervals, then checks if the true value of default (p1 for first class, p2 for second) is contained in the interval
### if it is, the joint probability of getting these estimates is added to a running tally to get the coverage probability.  If not, the coverage prob. is unaffected
### Next, the length of the interval is calclated and multiplied by the joint probability.  A running tally is kept to compute the expected length.

     ci <- ac(X1,N1,n1,alpha2)
	if ((p1>=ci[1]) && (p1<=ci[2])) {ac_cp1<-ac_cp1 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	ac_el1<-ac_el1+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)
     ci <- ac(X2,N2,n2,alpha2)
	if ((p2>=ci[1]) && (p2<=ci[2])) {ac_cp2<-ac_cp2 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	ac_el2<-ac_el2+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)

     #ci <- wald(X1,N1,n1,alpha2)
	#if ((p1>=ci[1]) && (p1<=ci[2])) {w_cp1<-w_cp1 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	#w_el1<-w_el1+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)
     #ci <- wald(X2,N2,n2,alpha2)
	#if ((p2>=ci[1]) && (p2<=ci[2])) {w_cp2<-w_cp2 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	#w_el2<-w_el2+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)

     ci <- score(X1,N1,n1,alpha2)
	if ((p1>=ci[1]) && (p1<=ci[2])) {s_cp1<-s_cp1 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	s_el1<-s_el1+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)
     ci <- score(X2,N2,n2,alpha2)
	if ((p2>=ci[1]) && (p2<=ci[2])) {s_cp2<-s_cp2 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	s_el2<-s_el2+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)

     ci <- jeff(X1,N1,n1,alpha2)
	if ((p1>=ci[1]) && (p1<=ci[2])) {j_cp1<-j_cp1 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	j_el1<-j_el1+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)
     ci <- jeff(X2,N2,n2,alpha2)
	if ((p2>=ci[1]) && (p2<=ci[2])) {j_cp2<-j_cp2 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	j_el2<-j_el2+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)

     ci <- zl(X1,N1,n1,alpha2)
	if ((p1>=ci[1]) && (p1<=ci[2])) {zl_cp1<-zl_cp1 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	zl_el1<-zl_el1+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)
     ci <- zl(X2,N2,n2,alpha2)
	if ((p2>=ci[1]) && (p2<=ci[2])) {zl_cp2<-zl_cp2 + dbinom(x1,n1,p1)*dbinom(x2,n2,p2)}
	zl_el2<-zl_el2+(ci[2]-ci[1])*dbinom(x1,n1,p1)*dbinom(x2,n2,p2)
		
  }
 }

### Once the loop completes, we have coverage probabilities and exptected lengths for each set of joint intervals.
### The following code returns a data frame of these values

cp<-c(ac_cp1,ac_cp2,w_cp1,w_cp2,s_cp1,s_cp2,j_cp1,j_cp2,zl_cp1,zl_cp2)
el<-c(ac_el1,ac_el2,w_el1,w_el2,s_el1,s_el2,j_el1,j_el2,zl_el1,zl_el2)
tab<-cbind(cp,el)
row.names(tab)<-c("ac1","ac2","w1","w2","s1","s2","j1","j2","zl1","zl2")
return(tab)

}


### ac calculates agresti-coull interval
### note that in the case of constraint violation, the original number of observations is still used to form the width of the interval

ac <- function(X,N,n,conf_lev){
z<-qnorm(1-conf_lev)
z2 <- z*z
xmod <- X + 0.5 * z2
nmod <- N + z2
p <- xmod/nmod
L <- p - z * sqrt(p * (1 - p)/n)
U <- p + z * sqrt(p * (1 - p)/n)
return(c(L,U))
}

### wald calculates the wald interval
### note that in the case of constraint violation, the original number of observations is still used to form the width of the interval

wald <- function(X,N,n,conf_lev){
z<-qnorm(1-conf_lev)
p <- X/N
L <- p - z * sqrt(p * (1 - p)/n)
U <- p + z * sqrt(p * (1 - p)/n)
return(c(L,U))
}

### score calculates the score interval
### note that in the case of constraint violation, the original number of observations is still used to form the width of the interval

score <- function(X,N,n,conf_lev){
z<-qnorm(1-conf_lev)
z2 <- z*z
p <- X/N
p1 <- p + 0.5 * z2/n
p2 <- z * sqrt((p * (1 - p) + 0.25 * z2/n)/n)
p3 <- 1 + z2/n
L <- (p1 - p2)/p3
U <- (p1 + p2)/p3
return(c(L,U))
}

### jeff calculates the jefferys interval

jeff <- function(X,N,n,conf_lev){
if (N != n) {X = X*n/N}	### Check if N is not equal to n; if so, there was a constraint violation, so we adjust X so that it is smaller and gives the same estimate when divided by the original number of observations, n.  Since Beta dist. can take positive, non-integer values this is ok.  
if (X==0) {L<-0} else   
{L <- qbeta(conf_lev, X+0.5, n-X+0.5)}
if (X==N)  {U=1} else
{U <- qbeta(1 - conf_lev, X+0.5, n-X+0.5)}
return(c(L,U))
}
    

### Function zl returns the lower and upper bounds of the Zhou-Li interval
### note that in the case of constraint violation, the original number of observations is still used to form the width of the interval

zl <- function(X,N,n,conf_lev){
if (X==0 || X==N) {X<-X+0.5; N<-N+1}
z<-qnorm(1-conf_lev)
phat<-X/N
qhat<-1-phat
a<- -1/6
b<- 1/6

if (phat==0.5) {gam<-2} else {gam<-(1-2*phat)/sqrt(phat*qhat)}   ### The value of gamma defined in paper by Zhou-Li is undefined for phat=0.5 but converges to 2 as phat-->0.5, so define to be 2.

ginv <- function(x)return(n^{1/2}*(a*gam)^{-1}*(sign(1+3*a*gam*(n^{-1/2}*x-n^{-1}*b*gam))*abs((1+3*a*gam*(n^{-1/2}*x-n^{-1}*b*gam)))^{1/3}-1))

L<-exp(log(phat/qhat)-(n*phat*qhat)^{-1/2}*ginv(z))/(1+exp(log(phat/qhat)-(n*phat*qhat)^{-1/2}*ginv(z)))
U<-exp(log(phat/qhat)-(n*phat*qhat)^{-1/2}*ginv(-z))/(1+exp(log(phat/qhat)-(n*phat*qhat)^{-1/2}*ginv(-z)))

return(c(L,U))
}
