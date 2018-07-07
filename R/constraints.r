heq <- function(x) {
h <- rep(NA,1)
h[1] <- sum(x)-1
h}

hin <- function(x) {
h <- rep(NA,1)
for (i in 1:length(x)) h[i] <- x[i]
h
}