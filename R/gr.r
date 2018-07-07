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