
########  Calculate remaining notional after j defaults ########
## input number defaults j,low attach point, al, and high attach point ah

rn <- function(j,al,ah) {
nl <- al*N_c/(1-R)
nh <- ah*N_c/(1-R)
if (j < ceiling(nl)) p <- (ah-al)*N_c*K
 else if (ceiling(nl) <= j & j < ceiling(nh)) p <- ah*N_c*K - j*(1-R)*K
 else if (j >= ceiling(nh)) p <- 0
return(p)
}