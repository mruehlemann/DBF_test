# DBF test
# --------

# Date:13/09/2012
# Author: Christopher Minas

"trace.matrix"<- function(x)
{
	return(sum(diag(x)))
}
"create.G" <- function(no.reps,D) 
{
	A <- -(1/2)*D^2
      ones <- matrix(data=1,ncol=no.reps,nrow=no.reps)
	centmat <- (diag(no.reps)-ones/no.reps)
      G <- centmat%*%A%*%centmat	### Hier
      return(G)
}

"perm.1st.moment" <- function(A,B,n)
{
	trace.matrix(A)*trace.matrix(B)/(n-1)
}

"perm.2nd.moment" <- function(A1,A2,A3,B1,B2,B3,n)
{
	p11 <- ((n-1)*A3*B3+(B1^2-B3)*(A1^2-A3)+2*(A2-A3)*(B2-B3)+4*A3*B3)/(n*(n-1))
	p22 <- (4*(n-3)*(2*A3-A2)*(2*B3-B2)+2*(2*A3-A1^2)*(2*B3-B1^2)*(n-3)+(2*A2+A1^2-6*A3)*(2*B2+B1^2-6*B3))/(n*(n-1)*(n-2)*(n-3))
	return(p11+p22)
}

"perm.3rd.moment" <- function(A1,A2,A3,A4,A5,A6,A7,A8,B1,B2,B3,B4,B5,B6,B7,B8,n)
{
	n2 <- n^2
	n3 <- n^3
	n4 <- n^4
	return((n2*(n+1)*(n2+15*n-4)*A5*B5 + 4*(n4-8*n3+19*n2-4*n-16)*A6*B6 + 24*(n2-n-4)*(A6*B8+B6*A8)+ 6*(n4-8*n3+21*n2-6*n-24)*A8*B8 + 12*(n4-n3-8*n2+36*n-48)*A7*B7 + 12*(n3-2*n2+9*n-12)*(A1*A3*B7 + A7*B1*B3) + 3*(n4 - 4*n3 - 2*n2+9*n-12)*A1*B1*A3*B3 + 24*((n3 - 3*n2 - 2*n+8)*(A7*B6 + A6*B7) + (n3 - 2*n2 - 3*n+12)*(A7*B8 + A8*B7)) + 12*(n2 - n + 4)*(A1*A3*B6 + B1*B3*A6) + 6*(2*n3 - 7*n2 - 3*n + 12)*(A1*A3*B8 + A8*B1*B3) - 2*n*(n-1)*(n2-n+4)*((2*A6+3*A8)*B5+(2*B6+3*B8)*A5) - 3*n*(n-1)*(n-1)*(n+4)*((A1*A3+4*A7)*B5+(B1*B3+4*B7)*A5) + 2*n*(n-1)*(n-2)*((A1^3 + 6*A1*A2 + 8*A4)*B5+(B1^3 + 6*B1*B2 + 8*B4)*A5) + (A1^3)*((n3-9*n2+23*n-14)*(B1^3)+6*(n-4)*B1*B2+8*B4) + 6*A1*A2*((n-4)*(B1^3)+(n3-9*n2+24*n-14)*B1*B2 + 4*(n-3)*B4) + 8*A4*((B1^3)+3*(n-3)*B1*B2+(n3-9*n2+26*n-22)*B4) - 16*((A1^3)*B6+A6*(B1^3)) - 6*(A1*A2*B6 + A6*B1*B2)*(2*n2-10*n+16) - 8*(A4*B6+A6*B4)*(3*n2-15*n+16)-((A1^3)*B8+A8*(B1^3))*(6*n2-30*n+24)-6*(A1*A2*B8+A8*B1*B2)*(4*n2-20*n+24) - 8*(A4*B8+A8*B4)*(3*n2-15*n+24) - (n-2)*(24*((A1^3)*B7+A7*(B1^3))+6*(A1*A2*B7+A7*B1*B2)*(2*n2-10*n+24)+8*(A4*B7+A7*B4)*(3*n2-15*n+24)+(3*n2-15*n+6)*((A1^3)*B1*B3+A1*A3*(B1^3))+6*(A1*A2*B1*B3+A1*A3*B1*B2)*(n2-5*n+6) + 48*(A4*B1*B3+A1*A3*B4)))/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)))
}

"asym.perm.info" <- function(A,B,n)
{
	Am2 <- A%*%A		### Hier
	Am3 <- A%*%Am2		### Hier 
	dA <- diag(A)
	Ap2 <- A^2
	Ap3 <- Ap2*A
	A1 <- sum(dA)
	A2 <- trace.matrix(Am2)
	A3 <- trace.matrix(Ap2)
	A4 <- trace.matrix(Am3)
	A5 <- trace.matrix(Ap3)
	A6 <- sum(Ap3)
	A7 <- dA%*%diag(Am2)
	A8 <- dA%*%A%*%dA
	Bm2 <- B%*%B		### Hier
	Bm3 <- B%*%Bm2		### Hier
	dB <- diag(B)
	Bp2 <- B^2
	Bp3 <- Bp2*B
	B1 <- sum(dB)
	B2 <- trace.matrix(Bm2)
	B3 <- trace.matrix(Bp2)
	B4 <- trace.matrix(Bm3)
	B5 <- trace.matrix(Bp3)
	B6 <- sum(B^3)
	B7 <- dB%*%diag(Bm2)
	B8 <- dB%*%B%*%dB
	mom1 <- perm.1st.moment(A, B, n)
    	mom2 <- perm.2nd.moment(A1,A2,A3,B1,B2,B3,n)
   	mom3 <- perm.3rd.moment(A1,A2,A3,A4,A5,A6,A7,A8,B1,B2,B3,B4,B5,B6,B7,B8,n)
   	mean.T <- mom1
      variance.T <- (mom2 - mom1^2)
      skewness.T <- (mom3 - 3 * variance.T * mom1 - mom1^3)/(variance.T^(3/2))
      tr.G <- trace.matrix(B)
	return(list(mean.T=mean.T,variance.T=variance.T,skewness.T=skewness.T,tr.G=tr.G))
}


"inv.f.fn" <- function(asym.info,f)
{
	if((f==-Inf)||(f==Inf)){
		ans <- ((asym.info$tr.G-asym.info$mean.T)-asym.info$mean.T/f)/(sqrt(asym.info$variance.T)*(1/f+1))
	}else{
		ans <-((asym.info$tr.G-asym.info$mean.T)*f-asym.info$mean.T)/(sqrt(asym.info$variance.T)*(1+f))
	}
	return(ans)
}
"fl" <- function(asym.info)
{
	a <- asym.info$tr.G
	m <- asym.info$mean.T
	s <- sqrt(asym.info$variance.T)
	g <- asym.info$skewness.T
	return(1/((g*a)/(g*m - 2*s)-1))
}

"cdf.F" <- function(asym.info,a)
{
	alpha <- fl(asym.info)
	if(asym.info$skewness.T >= 0){
		a1 <- pgamma(((asym.info$tr.G-asym.info$mean.T)/sqrt(asym.info$variance.T))+2/asym.info$skewness.T, shape = (4/(asym.info$skewness.T^2)), scale = (asym.info$skewness.T/2))
		if(a >= alpha){
			ans <- 1 + pgamma(inv.f.fn(asym.info,a)+2/asym.info$skewness.T, shape = (4/(asym.info$skewness.T^2)), scale = (asym.info$skewness.T/2)) - a1
		}
		if (a <= -1){
			ans <- pgamma(inv.f.fn(asym.info,a)+2/asym.info$skewness.T, shape = (4/(asym.info$skewness.T^2)), scale = (asym.info$skewness.T/2)) - a1
		}
		if( (-1 < a)&& (a <alpha)){
			ans <- 1 + pgamma(inv.f.fn(asym.info,alpha)+2/asym.info$skewness.T, shape = (4/(asym.info$skewness.T^2)), scale = (asym.info$skewness.T/2)) - a1
		}
	}
	if ((asym.info$skewness.T < 0)&&(alpha < -1)){
		a1 <- pgamma((2/abs(asym.info$skewness.T))-((asym.info$tr.G-asym.info$mean.T)/sqrt(asym.info$variance.T)), shape = (4/(asym.info$skewness.T^2)),scale = (abs(asym.info$skewness.T)/2),lower.tail=FALSE)
		if(a <= alpha){
			ans <- pgamma((2/abs(asym.info$skewness.T))-inv.f.fn(asym.info,a), shape = (4/(asym.info$skewness.T^2)), scale = (abs(asym.info$skewness.T)/2),lower.tail=FALSE) - a1
		}
		if (a > -1){
			ans <- 1 + pgamma((2/abs(asym.info$skewness.T))-inv.f.fn(asym.info,a), shape = (4/(asym.info$skewness.T^2)), scale = (abs(asym.info$skewness.T)/2),lower.tail=FALSE) - a1
		}
		if( (alpha < a)&& (a <= -1)){
			ans <- 1 + pgamma((2/abs(asym.info$skewness.T))-inv.f.fn(asym.info,-1), shape = (4/(asym.info$skewness.T^2)), scale = (abs(asym.info$skewness.T)/2),lower.tail=FALSE) - a1
		}
	}
	if ((asym.info$skewness.T < 0)&&(alpha > -1)){
		a1 <- pgamma((2/abs(asym.info$skewness.T))-((asym.info$tr.G-asym.info$mean.T)/sqrt(asym.info$variance.T)), shape = (4/(asym.info$skewness.T^2)),scale = (abs(asym.info$skewness.T)/2),lower.tail=FALSE)
		if(a < -1){
			ans <- 0
		}
		if((-1 < a)&&(a<= alpha)){
			ans <- pgamma((2/abs(asym.info$skewness.T))-inv.f.fn(asym.info,a), shape = (4/(asym.info$skewness.T^2)),scale = (abs(asym.info$skewness.T)/2),lower.tail=FALSE)
		}
		if(a > alpha){
			ans <- 1
		}
	}	
	return(ans)
}

"pearson.type.3.dgamma" <- function(g,y)
{
	if(g >= 0){
		ans <- dgamma(y+2/g, shape=4/(g^2), scale = g/2, log = FALSE)
	}else{
		ans <- dgamma(abs(2/g) -y, shape=4/(g^2), scale = abs(g)/2, log = FALSE)
	}
	return(ans)
}
"pdf.F" <- function(asym.info,ys)
{
	a <- asym.info$tr.G
	m <- asym.info$mean.T
	s <- sqrt(asym.info$variance.T)
	g <- asym.info$skewness.T
	if(ys==-1){
		vals <- 0
	}else{
		vals <- a/(s*(1+ys)^2)*pearson.type.3.dgamma(g,inv.f.fn(asym.info,ys))
	}
	return(vals)
}


"DBF.stat.G" <- function(G,H)
{
	return((1/((trace.matrix(G)/trace.matrix(H%*%G))-1)))
}

"make.H" <- function(N,group.labels)
{
	unique.grp.labels <- unique(group.labels)
	no.grps <- length(unique.grp.labels)
	ind.grps <- lapply(1:no.grps, function(i) which(group.labels==unique.grp.labels[i]))
	L <- matrix(0,nrow=N,ncol=length(unique.grp.labels))
	for (i in 1:no.grps){
		L[ind.grps[[i]],i] <- 1
	}	
	return(L%*%solve(t(L)%*%L)%*%t(L))
}


### DBF.test ist die Hauptfunktion die aufgerufen wird
"DBF.test" <- function(dmat,group.labels,n)
{
	gmat <- create.G(n,dmat)
	H <- scale(make.H(n,group.labels),center=TRUE,scale=FALSE)
	dbf0 <- DBF.stat.G(gmat,H)
	p.dbf0 <- 1- cdf.F(asym.perm.info(H,gmat,n),dbf0)
	return(c(dbf.statistic=dbf0,dbf.p.value=p.dbf0))
}




############## Example:
# --------

"rmnorm" <- function (n = 1, mean = rep(0, d), varcov) # returns multivariate random normal observation
{
    d <- if (is.matrix(varcov))
        ncol(varcov)
    else 1
    z <- matrix(rnorm(n * d), n, d) %*% chol(varcov)
    y <- t(mean + t(z))
    return(y)
}

"rwishart" <- function(df, p) # returns random Wishart matrix
{
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, df:(df-p+1)))
    if(p > 1)
    {
        pseq <- 1:(p-1)
        Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    return(crossprod(Z))
}

library(MASS) 
library(graphics)
library(stats)
library(Matrix)
library(combinat)
library(ecodist)

set.seed(101)

# create multivariate normal dataset of size 99 x 10:

n <- 99
p <- 10
varcov <- rwishart(df=p,p=p)
mean1 <- runif(p,-6,6)
X <- rmnorm(n, mean=mean1, varcov)	

# set 3 equal-sized groups:

n1 <- n2 <- n3 <- n/3
group.labels <- c(rep(1,n1),rep(2,n2),rep(3,n3))	

# define Mahalanobis distance matrix:

dmat <- as.matrix(distance(X,"mahalanobis"))

# apply DBF test of equality between groups:
DBF.test(dmat,group.labels,n)

#dbf.statistic   dbf.p.value
#   0.01487275    0.65814965


## plot approximate CDF of DBF:

asym.info <- asym.perm.info(A=scale(make.H(n,group.labels),center=TRUE,scale=FALSE),B=create.G(n,dmat),n)
xs <- seq(fl(asym.info),0.2,len=1000)	
cdf.dbf <- sapply(1:1000, function(i) cdf.F(asym.info,xs[i]))
plot(xs,cdf.dbf,type="l")

## plot approximate PDF of DBF:

xs <- seq(fl(asym.info),0.2,len=1000)	
pdf.dbf <- sapply(1:1000, function(i) pdf.F(asym.info,xs[i]))
plot(xs,pdf.dbf,type="l")

## comparison of approximate PDF with sampling distribution using 10^4 Monte Carlo permutations:

gmat <- create.G(n,dmat)
H <- scale(make.H(n,group.labels),center=TRUE,scale=FALSE)
dbf.perms <- sapply(1:10000, function(i) {perms <- sample(1:n,replace=FALSE)
							DBF.stat.G(G=gmat[perms,perms],H)})
plot(xs,pdf.dbf,type="l")
hist(dbf.perms,prob=TRUE,add=TRUE)
