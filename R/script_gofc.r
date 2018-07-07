a <- .4
n <- 200
len <- 10000
cloppercp<-rep(0,len)
clopperel<-rep(0,len)
wcp<-rep(0,len)
wel<-rep(0,len)
w2cp<-rep(0,len)
w2el<-rep(0,len)
accp<-rep(0,len)
acel<-rep(0,len)
ac2cp<-rep(0,len)
ac2el<-rep(0,len)
scp<-rep(0,len)
sel<-rep(0,len)
s2cp<-rep(0,len)
s2el<-rep(0,len)
mscp<-rep(0,len)
msel<-rep(0,len)
zlcp<-rep(0,len)
zlel<-rep(0,len)
zl2cp<-rep(0,len)
zl2el<-rep(0,len)
jeffcp<-rep(0,len)
jeffel<-rep(0,len)
mjcp<-rep(0,len)
mjel<-rep(0,len)
jeffcp<-rep(0,len)
jeffel<-rep(0,len)
exactcp<-rep(0,len)
exactel<-rep(0,len)
clopper_temp <- clopper_gofc(n, .5)
clopperL <- clopper_temp[1:(n+1)]
clopperU <- clopper_temp[(n+2):(2*n+2)]
w_temp <- wald_gofc(n, .5)
wL <- w_temp[1:(n+1)]
wU <- w_temp[(n+2):(2*n+2)]
wL2 <- w_temp[(2*n+3):(3*n+3)]
wU2 <- w_temp[(3*n+4):(4*n+4)]
ac_temp <- ac_gofc(n, .5)
acL <- ac_temp[1:(n+1)]
acU <- ac_temp[(n+2):(2*n+2)]
acL2 <- ac_temp[(2*n+3):(3*n+3)]
acU2 <- ac_temp[(3*n+4):(4*n+4)]
score_temp <- score_gofc(n, .5)
scoreL <- score_temp[1:(n+1)]
scoreU <- score_temp[(n+2):(2*n+2)]
scoreL2 <- score_temp[(2*n+3):(3*n+3)]
scoreU2 <- score_temp[(3*n+4):(4*n+4)]
jeff_temp <- jeff_gofc(n, .5)
jeffL <- jeff_temp[1:(n+1)]
jeffU <- jeff_temp[(n+2):(2*n+2)]
ms_temp <- mod_score_gofc(n, .5)
msL <- ms_temp[1:(n+1)]
msU <- ms_temp[(n+2):(2*n+2)]
mj_temp <- mod_jeff_gofc(n, .5)
mjL <- mj_temp[1:(n+1)]
mjU <- mj_temp[(n+2):(2*n+2)]
zl_temp <- zl_gofc(n, .5)
zlL <- zl_temp[1:(n+1)]
zlU <- zl_temp[(n+2):(2*n+2)]
zlL2 <- zl_temp[(2*n+3):(3*n+3)]
zlU2 <- zl_temp[(3*n+4):(4*n+4)]
exact_temp <- exact_gofc(n, .5)
exactL <- exact_temp[1:(n+1)]
exactU <- exact_temp[(n+2):(2*n+2)]

points <- seq(.00004,.5,length=len)

for (i in 1:len) cloppercp[i] <- cp(n,points[i],clopperL,clopperU,p_gofc_200_4[[i]])
for (i in 1:len) clopperel[i] <- el(n,points[i],clopperL,clopperU,p_gofc_200_4[[i]])
for (i in 1:len) wcp[i] <- cp(n,points[i],wL,wU,p_gofc_200_4[[i]])
for (i in 1:len) wel[i] <- el(n,points[i],wL,wU,p_gofc_200_4[[i]])
for (i in 1:len) w2cp[i] <- cp(n,points[i],wL2,wU2,p_gofc_200_4[[i]])
for (i in 1:len) w2el[i] <- el(n,points[i],wL2,wU2,p_gofc_200_4[[i]])
for (i in 1:len) accp[i] <- cp(n,points[i],acL,acU,p_gofc_200_4[[i]])
for (i in 1:len) acel[i] <- el(n,points[i],acL,acU,p_gofc_200_4[[i]])
for (i in 1:len) ac2cp[i] <- cp(n,points[i],acL2,acU2,p_gofc_200_4[[i]])
for (i in 1:len) ac2el[i] <- el(n,points[i],acL2,acU2,p_gofc_200_4[[i]])
for (i in 1:len) scp[i] <- cp(n,points[i],scoreL,scoreU,p_gofc_200_4[[i]])
for (i in 1:len) sel[i] <- el(n,points[i],scoreL,scoreU,p_gofc_200_4[[i]])
for (i in 1:len) s2cp[i] <- cp(n,points[i],scoreL2,scoreU2,p_gofc_200_4[[i]])
for (i in 1:len) s2el[i] <- el(n,points[i],scoreL2,scoreU2,p_gofc_200_4[[i]])
for (i in 1:len) mscp[i] <- cp(n,points[i],msL,msU,p_gofc_200_4[[i]])
for (i in 1:len) msel[i] <- el(n,points[i],msL,msU,p_gofc_200_4[[i]])
for (i in 1:len) mjcp[i] <- cp(n,points[i],mjL,mjU,p_gofc_200_4[[i]])
for (i in 1:len) mjel[i] <- el(n,points[i],mjL,mjU,p_gofc_200_4[[i]])
for (i in 1:len) zlcp[i] <- cp(n,points[i],zlL,zlU,p_gofc_200_4[[i]])
for (i in 1:len) zlel[i] <- el(n,points[i],zlL,zlU,p_gofc_200_4[[i]])
for (i in 1:len) zl2cp[i] <- cp(n,points[i],zlL2,zlU2,p_gofc_200_4[[i]])
for (i in 1:len) zl2el[i] <- el(n,points[i],zlL2,zlU2,p_gofc_200_4[[i]])
for (i in 1:len) jeffcp[i] <- cp(n,points[i],jeffL,jeffU,p_gofc_200_4[[i]])
for (i in 1:len) jeffel[i] <- el(n,points[i],jeffL,jeffU,p_gofc_200_4[[i]])
for (i in 1:len) exactcp[i] <- cp(n,points[i],exactL,exactU,p_gofc_200_4[[i]])
for (i in 1:len) exactel[i] <- el(n,points[i],exactL,exactU,p_gofc_200_4[[i]])

n200cloppercp4 <- cloppercp
n200clopperel4 <- clopperel
n200wcp4 <- wcp
n200wel4 <- wel
n200w2cp4 <- w2cp
n200w2el4 <- w2el
n200accp4 <- accp
n200acel4 <- acel
n200ac2cp4 <- ac2cp
n200ac2el4 <- ac2el
n200scp4 <- scp
n200sel4 <- sel
n200s2cp4 <- s2cp
n200s2el4 <- s2el
n200mscp4 <- mscp
n200msel4 <- msel
n200zlcp4 <- zlcp
n200zlel4 <- zlel
n200zl2cp4 <- zl2cp
n200zl2el4 <- zl2el
n200jeffcp4 <- jeffcp
n200jeffel4 <- jeffel
n200mjcp4 <- mjcp
n200mjel4 <- mjel
n200jeffcp4 <- jeffcp
n200jeffel4 <-  jeffel
n200exactcp4 <- exactcp
n200exactel4 <- exactel


plot(points,n100exactcp1,xlab="probability",ylab="CP",ylim=c(.94,1),type="l",lty=1,col="blue")
lines(points,n100exactcp2,type="l",lty=2,col="red")
lines(points,n100exactcp3,type="l",lty=3,col="black")
lines(points,n100exactcp4,type="l",lty=4,col="cyan")
lines(points,exactcp,type="l",lty=1,col="black")
abline(h=.95)
title("Coverage Prob., Exact, n=100")
#title(bquote(bold(paste("Coverage Probability, n=",.(n)))))
#legend(x=.2,y=.925,c("AC2","Score","ZL2","Jefferys","Exact"),cex=.75,col=c("blue","red","cyan","green","black"),lty=c(1,2,3,4,1))
legend(x=.4,y=.99,c("a=0.1","a=0.2","a=0.3","a=0.4"),cex=.75,col=c("blue","red","black","cyan"),lty=c(1,2,3,4))
legend(x=.4,y=.99,c("n=50","n=100","n=200"),cex=.75,col=c("blue","red","black"),lty=c(1,2,3))
legend(x=.4,y=.6,c("ZL","ZL2"),cex=.75,col=c("blue","red"),lty=c(1,2))
mtext("correlation = 0.1,0.2,0.3,0.4")
mtext("correlation = 0.2")

plot(points,n100exactel1,xlab="probability",ylab="EL",ylim=c(0.04,.6),type="l",lty=1,col="blue")
lines(points,n100exactel2,type="l",lty=2,col="red")
lines(points,n100exactel3,type="l",lty=3,col="black")
lines(points,n100exactel4,type="l",lty=4,col="cyan")
lines(points,exactel,type="l",lty=1,col="black")
title("Expected Length, Exact, n=100")
mtext("correlation = 0.2")
#legend(x=.2,y=.1,c("AC2","Score","ZL2","Jefferys","Exact"),cex=.75,col=c("blue","red","cyan","green","black"),lty=c(1,2,3,4,1))
legend(x=.3,y=.1,c("n=50","n=100","n=200"),cex=.75,col=c("blue","red","black"),lty=c(1,2,3))
legend(x=.4,y=.15,c("a=0.1","a=0.2","a=0.3","a=0.4"),cex=.75,col=c("blue","red","black","cyan"),lty=c(1,2,3,4))
legend(x=.3,y=.075,c("Wald","Wald2"),cex=.75,col=c("blue","red"),lty=c(1,2))

#### some constants for faster testing
p_gofc_50_1 <- p_gofc_50_2 <- p_gofc_50_3 <- p_gofc_50_4 <- 
p_gofc_100_1 <-p_gofc_100_2 <-p_gofc_100_3 <-p_gofc_100_4 <-
p_gofc_200_1 <-p_gofc_200_2 <-p_gofc_200_3 <-p_gofc_200_4 <- 
list()

a<-.1
for (i in 1:len) p_gofc_200_1[[i]] <- gofc_prob_vec(200,points[i])
for (i in 1:len) p_gofc_100_1[[i]] <- gofc_prob_vec(100,points[i])
for (i in 1:len) p_gofc_50_1[[i]] <- gofc_prob_vec(50,points[i])
a<-.2
for (i in 1:len) p_gofc_200_2[[i]] <- gofc_prob_vec(200,points[i])
for (i in 1:len) p_gofc_100_2[[i]] <- gofc_prob_vec(100,points[i])
for (i in 1:len) p_gofc_50_2[[i]] <- gofc_prob_vec(50,points[i])
a<-.3
for (i in 1:len) p_gofc_200_3[[i]] <- gofc_prob_vec(200,points[i])
for (i in 1:len) p_gofc_100_3[[i]] <- gofc_prob_vec(100,points[i])
for (i in 1:len) p_gofc_50_3[[i]] <- gofc_prob_vec(50,points[i])
a<-.4
for (i in 1:len) p_gofc_200_4[[i]] <- gofc_prob_vec(200,points[i])
for (i in 1:len) p_gofc_100_4[[i]] <- gofc_prob_vec(100,points[i])
for (i in 1:len) p_gofc_50_4[[i]] <- gofc_prob_vec(50,points[i])



