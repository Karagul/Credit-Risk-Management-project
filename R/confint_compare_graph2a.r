### function confint_compare_graph takes input of: probability, p, n1, n2, conf.level,
### fix, which can equal 1 or 2 based on whether we want to fix p1 or p2 as p, and
### step, which is the step size that the not fixed probability increases by as we compute
### coverage prob and expected length.
### The function returns graphs of CP and EL for both intervals

confint_compare_graph2<-function(p,n1,n2,conf.level, fix=2, step=.001) {

par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
op<-par(mfrow = c(2, 2), oma = c(0, 0, 2, 0)) 

i=1;
if (fix==1) {
	p1<-p;
	j<-seq(p1+step, .5, step);
	tab<-matrix(nrow=20,ncol=length(j))
	for (p2 in j) {
		tab[,i]<-as.vector(confint_compare2(p1,n1,p2,n2,conf.level))
		i<-i+1;
		}
	}
else {
	p2<-p;
	j<-seq(step,p2-step,step);
	tab<-matrix(nrow=20,ncol=length(j))
	for (p1 in j) {
		tab[,i]<-as.vector(confint_compare2(p1,n1,p2,n2,conf.level))
		i<-i+1;
		}
	}

plot(j,tab[1,],xlab="prob.",ylab="CP",ylim=c(.9,1),
#ylim=c(min(tab[1,],tab[5,],tab[7,],tab[9,]),max(tab[1,],tab[3,],tab[5,],tab[7,],tab[9,])),
type="l",lty=2,col="red",main="Interval 1 Coverage Probability");
abline(h=.95)
lines(j,tab[5,],lty=3,col="blue");
lines(j,tab[7,],lty=4,col="green");
lines(j,tab[9,],lty=5,col="cyan");


plot(j,tab[2,],xlab="prob.",ylab="CP",ylim=c(.9,1),
#ylim=c(min(tab[2,],tab[6,],tab[8,],tab[10,]),max(tab[2,],tab[4,],tab[6,],tab[8,],tab[10,])),
type="l",lty=2,col="red",main="Interval 2 Coverage Probability")
abline(h=.95)
lines(j,tab[6,],lty=3,col="blue")
lines(j,tab[8,],lty=4,col="green")
lines(j,tab[10,],lty=5,col="cyan")


plot(j,tab[11,],xlab="prob.",ylab="EL",
ylim=c(min(tab[11,],tab[15,],tab[17,],tab[19,]),max(tab[11,],tab[13,],tab[15,],tab[17,],tab[19,])),
type="l",lty=2,col="red",main="Interval 1 Expected Length")
lines(j,tab[15,],lty=3,col="blue")
lines(j,tab[17,],lty=4,col="green")
lines(j,tab[19,],lty=5,col="cyan")


plot(j,tab[12,],xlab="prob.",ylab="EL",
ylim=c(min(tab[12,],tab[16,],tab[18,],tab[20,]),max(tab[12,],tab[14,],tab[16,],tab[18,],tab[20,])),
type="l",lty=2,col="red",main="Interval 2 Expected Length")
lines(j,tab[16,],lty=3,col="blue")
lines(j,tab[18,],lty=4,col="green")
lines(j,tab[20,],lty=5,col="cyan")


par(op)

if (fix==1) {title(bquote(paste("Strat2: p2=prob.; ",p1==.(p),"; ",n1==.(n1),"; ",n2==.(n2),"; ",step==.(step))),outer=T)}
else {title(bquote(paste("Strat2: p1=prob.; ",p2==.(p),"; ",n1==.(n1),"; ",n2==.(n2),"; ",step==.(step))),outer=T)}

op <- par(usr=c(0,1,0,1), xpd=NA)
legend(x=-.2,y=1.45, c(".95","AC","Score","Jefferys","ZL"), cex=0.6,
col=c("black","red","blue","green","cyan"), lty=c(1,2,3,4,5))


}

	
	


		
