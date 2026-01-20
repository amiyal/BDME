rm(list=ls(all=TRUE))
library(maxLik)
library(VGAM)
## Data -----------------------2
x=c(29, 25, 50, 15, 13, 27, 15, 18, 7, 7, 8, 19, 12, 18, 5, 21, 
15, 86, 21, 15, 14, 39, 15, 14, 70, 44, 6, 23, 58, 19, 50, 23, 
11, 6, 34, 18, 28, 34, 12, 37, 4, 60, 20, 23, 40, 65, 19, 31)


var(x)/mean(x) 
#ID--------
n=length(x)
xbar <- mean(x)


############################################################
######################## BDME Model ########################
############################################################
loglik.bdme <- function(theta) {
  p <- theta[1]
  if (p <= 0 || p >= 1) return(-Inf)
  2 * n * log(1 - p) +
    (2 * n * xbar - n) * log(p) +
    sum(log(x + p * (1 - p))) -
    (n * xbar + 2 * n) * log(1 - p * (1 - p))
}

ml.bdme <- maxLik(loglik.bdme, start=c(.9), method="NR")
p.bdme  <- ml.bdme$estimate
se.bdme <- sqrt(diag(vcov(ml.bdme)))

Fx.bdme <- function(x){
  1 - (p.bdme^(2*x+1))*(1+(1-p.bdme)*(x-p.bdme^2))/((1-p.bdme*(1-p.bdme))^(x+2))
}

ks.bdme  <- ks.test(x,"Fx.bdme")
AIC.bdme <- -2*logLik(ml.bdme)+2
BIC.bdme <- -2*logLik(ml.bdme)+log(n)

############################################################
######################## Poisson ###########################
############################################################
loglik.pois <- function(theta){
  lambda <- theta[1]
  if(lambda<=0) return(-Inf)
  sum(dpois(x, lambda, log=TRUE))
}

ml.pois <- maxLik(loglik.pois, start=c(xbar))
lam.hat <- ml.pois$estimate
se.pois <- sqrt(diag(vcov(ml.pois)))

Fx.pois <- function(x) ppois(x, lam.hat)

ks.pois  <- ks.test(x,"Fx.pois")
AIC.pois <- -2*logLik(ml.pois)+2
BIC.pois <- -2*logLik(ml.pois)+log(n)

############################################################
######################## Binomial ##########################
############################################################
m <- max(x)

loglik.bin <- function(theta){
  p <- theta[1]
  if(p<=0 || p>=1) return(-Inf)
  sum(dbinom(x, size=m, prob=p, log=TRUE))
}

ml.bin <- maxLik(loglik.bin, start=c(.5))
p.bin  <- ml.bin$estimate
se.bin <- sqrt(diag(vcov(ml.bin)))

Fx.bin <- function(x) pbinom(x, size=m, prob=p.bin)

ks.bin  <- ks.test(x,"Fx.bin")
AIC.bin <- -2*logLik(ml.bin)+2
BIC.bin <- -2*logLik(ml.bin)+log(n)

############################################################
################ Shifted Geometric #########################
############################################################
loglik.sgeom <- function(theta){
  p <- theta[1]
  if(p<=0 || p>=1) return(-Inf)
  sum(log(p) + (x-1)*log(1-p))
}

ml.sgeom <- maxLik(loglik.sgeom, start=c(.5))
p.sg  <- ml.sgeom$estimate
se.sg <- sqrt(diag(vcov(ml.sgeom)))

Fx.sgeom <- function(x){
  1 - (1 - p.sg)^x
}

ks.sgeom  <- ks.test(x,"Fx.sgeom")
AIC.sgeom <- -2*logLik(ml.sgeom)+2
BIC.sgeom <- -2*logLik(ml.sgeom)+log(n)

############################################################
################ Negative Binomial #########################
############################################################
loglik.nb <- function(theta){
  r <- theta[1]
  p <- theta[2]
  if(r<=0 || p<=0 || p>=1) return(-Inf)
  sum(dnbinom(x, size=r, prob=p, log=TRUE))
}

ml.nb <- maxLik(loglik.nb, start=c(1,0.5))
r.nb  <- ml.nb$estimate[1]
p.nb  <- ml.nb$estimate[2]
se.nb <- sqrt(diag(vcov(ml.nb)))

Fx.nb <- function(x){
  pnbinom(x, size=r.nb, prob=p.nb)
}

ks.nb  <- ks.test(x,"Fx.nb")
AIC.nb <- -2*logLik(ml.nb)+2*2
BIC.nb <- -2*logLik(ml.nb)+2*log(n)

############################################################
#################### Discrete Perks ########################
############################################################
loglik.dp <- function(theta){
  alpha <- theta[1]
  beta  <- theta[2]
  if(alpha<=0 || beta<=0) return(-Inf)

  sum(log(alpha) + log(1+alpha) + log(exp(beta)-1) + beta*x -
        log(1 + alpha*exp(beta*x)) -
        log(1 + alpha*exp(beta*(x+1))))
}

ml.dp <- maxLik(loglik.dp, start=c(0.44,0.08), method="NR")

alpha.dp <- ml.dp$estimate[1]
beta.dp  <- ml.dp$estimate[2]
se.dp    <- sqrt(diag(vcov(ml.dp)))

Fx.dp <- function(x){
  alpha.dp*(exp(beta.dp*(x+1))-1)/(1+alpha.dp*exp(beta.dp*(x+1)))
}

ks.dp  <- ks.test(x,"Fx.dp")
AIC.dp <- -2*logLik(ml.dp)+2*2
BIC.dp <- -2*logLik(ml.dp)+2*log(n)

############################################################
#################### Discrete Weibull ######################
############################################################
loglik.dw <- function(theta){
  beta <- theta[1]
  p    <- theta[2]
  if(beta<=0 || p<=0 || p>=1) return(-Inf)

  sum(log((1-p)^(x^beta) - (1-p)^((x+1)^beta)))
}

ml.dw <- maxLik(loglik.dw, start=c(1,0.5), method="NR")

beta.dw <- ml.dw$estimate[1]
p.dw    <- ml.dw$estimate[2]
se.dw   <- sqrt(diag(vcov(ml.dw)))

Fx.dw <- function(x){
  1 - (1 - p.dw)^((x+1)^beta.dw)
}

ks.dw  <- ks.test(x,"Fx.dw")
AIC.dw <- -2*logLik(ml.dw)+2*2
BIC.dw <- -2*logLik(ml.dw)+2*log(n)

############################################################
################ Generalized Poisson #######################
############################################################
loglik.gp <- function(theta){
  lambda1 <- theta[1]
  lambda2 <- theta[2]
  if(lambda1<=0 || abs(lambda2)>=1) return(-Inf)

  sum(dgenpois(x, lambda1, lambda2, log=TRUE))
}

# helper: pmf (Consul & Jain form)
dgenpois <- function(x, lambda1, lambda2, log=FALSE){
  val <- lambda1 * (lambda1 + lambda2*x)^(x-1) * 
         exp(-(lambda1 + lambda2*x)) / factorial(x)
  if(log) return(log(val)) else return(val)
}

# helper: cdf (numerical)
pgenpois <- function(q, lambda1, lambda2){
  sapply(q, function(k) sum(dgenpois(0:k, lambda1, lambda2)))
}

ml.gp <- maxLik(loglik.gp, start=c(mean(x), 0.1), method="NR")

lambda1.gp <- ml.gp$estimate[1]
lambda2.gp <- ml.gp$estimate[2]
se.gp      <- sqrt(diag(vcov(ml.gp)))

Fx.gp <- function(x){
  pgenpois(x, lambda1.gp, lambda2.gp)
}

ks.gp  <- ks.test(x,"Fx.gp")
AIC.gp <- -2*logLik(ml.gp)+2*2
BIC.gp <- -2*logLik(ml.gp)+2*log(n)

############################################################
################### Discrete Logistic ######################
############################################################

# CDF (CORRECT)
Fdlog <- function(x, p, mu){
  1/(1 + p^(x - mu + 1))
}

# PMF (difference of CDFs – CORRECT)
ddlog <- function(x, p, mu){
  Fdlog(x, p, mu) - Fdlog(x-1, p, mu)
}

# Log-likelihood
loglik.dlog <- function(theta){
  p  <- theta[1]
  mu <- theta[2]
  if(p<=0 || p>=1) return(-Inf)

  sum(log(ddlog(x, p, mu)))
}

# MLE
ml.dlog <- maxLik(loglik.dlog, start=c(0.9, 23), method="NR")

p.dlog  <- ml.dlog$estimate[1]
mu.dlog <- ml.dlog$estimate[2]
se.dlog <- sqrt(diag(vcov(ml.dlog)))

# CDF for KS
Fx.dlog <- function(x){
  Fdlog(x, p.dlog, mu.dlog)
}

# Goodness-of-fit
ks.dlog  <- ks.test(x, "Fx.dlog")
AIC.dlog <- -2*logLik(ml.dlog) + 2*2
BIC.dlog <- -2*logLik(ml.dlog) + 2*log(n)

############################################################
#################### Discrete Burr #########################
############################################################

# CDF
Fburr <- function(x, theta, alpha){
  1 - theta^(log(1 + (x+1)^alpha))
}

# PMF = difference of CDFs
dburr <- function(x, theta, alpha){
  Fburr(x, theta, alpha) - Fburr(x-1, theta, alpha)
}

# Log-likelihood
loglik.dburr <- function(par){
  theta <- par[1]
  alpha <- par[2]
  if(theta<=0 || theta>=1 || alpha<=0) return(-Inf)

  sum(log(dburr(x, theta, alpha)))
}

# MLE
ml.dburr <- maxLik(loglik.dburr, start=c(0.7,5), method="NM")

theta.db <- ml.dburr$estimate[1]
alpha.db <- ml.dburr$estimate[2]
se.dburr <- sqrt(diag(vcov(ml.dburr)))

# CDF for KS
Fx.dburr <- function(x){
  Fburr(x, theta.db, alpha.db)
}

# Goodness-of-fit
ks.dburr  <- ks.test(x, "Fx.dburr")
AIC.dburr <- -2*logLik(ml.dburr) + 2*2
BIC.dburr <- -2*logLik(ml.dburr) + 2*log(n)

############################################################
#################### Discrete Gamma ########################
############################################################

# Survival function
Sgamma <- function(x, shape, scale){
  1 - pgamma(x, shape=shape, scale=scale)
}

# PMF
dgamma.disc <- function(x, shape, scale){
  Sgamma(x, shape, scale) - Sgamma(x+1, shape, scale)
}

# Log-likelihood
loglik.dg <- function(par){
  shape <- par[1]
  scale <- par[2]
  if(shape<=0 || scale<=0) return(-Inf)

  sum(log(dgamma.disc(x, shape, scale)))
}

# MLE
ml.dg <- maxLik(loglik.dg, start=c(1, mean(x)), method="NR")

shape.dg <- ml.dg$estimate[1]
scale.dg <- ml.dg$estimate[2]
se.dg    <- sqrt(diag(vcov(ml.dg)))

# CDF for KS
Fx.dg <- function(x){
  1 - Sgamma(x+1, shape.dg, scale.dg)
}

# Goodness-of-fit
ks.dg  <- ks.test(x, "Fx.dg")
AIC.dg <- -2*logLik(ml.dg) + 2*2
BIC.dg <- -2*logLik(ml.dg) + 2*log(n)

############################################################
################# Discrete Lindley II by Hasan et al 2014###
############################################################

# PMF
dlin2 <- function(x, theta){
  p <- exp(-theta)
  (p^x/(1+theta)) * (theta*(1-2*p) + (1-p)*(1 + theta*x))
}

# CDF
Flin2 <- function(x, theta){
  p <- exp(-theta)
  1 - ((1 + theta + theta*x)/(1 + theta)) * p^x
}

# Log-likelihood
loglik.dlin2 <- function(theta){
  if(theta <= 0) return(-Inf)
  sum(log(dlin2(x, theta)))
}

# MLE
ml.dlin2 <- maxLik(loglik.dlin2, start=c(1), method="NR")

theta.dlin2 <- ml.dlin2$estimate
se.dlin2    <- sqrt(diag(vcov(ml.dlin2)))

# CDF for KS
Fx.dlin2 <- function(x){
  Flin2(x, theta.dlin2)
}

# Goodness-of-fit
ks.dlin2  <- ks.test(x, "Fx.dlin2")
AIC.dlin2 <- -2*logLik(ml.dlin2) + 2
BIC.dlin2 <- -2*logLik(ml.dlin2) + log(n)

############################################################
######## Exponentiated Discrete Lindley (EDLi) ##############
######## El-Morshedy et al. (2020) ##########################
############################################################

# Lambda function (as in paper)
Lambda <- function(x, a, b){
  (1 - a^x + ((1 + x)*a^x - 1)*log(a))^b
}

# CDF (Eq. 6)
Fedl <- function(x, a, b){
  Lambda(x+1, a, b) / (1 - log(a))^b
}

# PMF (Eq. 7)
dedl <- function(x, a, b){
  (Lambda(x+1, a, b) - Lambda(x, a, b)) / (1 - log(a))^b
}

# Log-likelihood
loglik.edl <- function(par){
  a <- par[1]
  b <- par[2]
  if(a <= 0 || a >= 1 || b <= 0) return(-Inf)
  sum(log(dedl(x, a, b)))
}

# MLE
ml.edl <- maxLik(loglik.edl, start=c(0.91, 1.34), method="NR")

a.edl  <- ml.edl$estimate[1]
b.edl  <- ml.edl$estimate[2]
se.edl <- sqrt(diag(vcov(ml.edl)))

# CDF for KS
Fx.edl <- function(x){
  Fedl(x, a.edl, b.edl)
}

# Goodness-of-fit
ks.edl  <- ks.test(x, "Fx.edl")
AIC.edl <- -2*logLik(ml.edl) + 2*2
BIC.edl <- -2*logLik(ml.edl) + 2*log(n)

############################################################
######## Exponentiated Geometric (EGD) #####################
######## Chakraborty & Gupta (2015) ########################
############################################################

# CDF (from paper)
Feg <- function(x, q, alpha){
  (1 - q^(x+1))^alpha
}

# PMF (Eq. 2)
deg <- function(x, q, alpha){
  (1 - q^(x+1))^alpha - (1 - q^x)^alpha
}

# Log-likelihood
loglik.eg <- function(par){
  q     <- par[1]
  alpha <- par[2]
  if(q<=0 || q>=1 || alpha<=0) return(-Inf)
  sum(log(deg(x, q, alpha)))
}

# MLE
ml.eg <- maxLik(loglik.eg, start=c(0.93, 2.68), method="NR")

q.eg     <- ml.eg$estimate[1]
alpha.eg <- ml.eg$estimate[2]
se.eg    <- sqrt(diag(vcov(ml.eg)))

# CDF for KS
Fx.eg <- function(x){
  Feg(x, q.eg, alpha.eg)
}

# Goodness-of-fit
ks.eg  <- ks.test(x, "Fx.eg")
AIC.eg <- -2*logLik(ml.eg) + 2*2
BIC.eg <- -2*logLik(ml.eg) + 2*log(n)

############################################################
######## Two Parameter Discrete Lindley (TDL) ##############
######## Hussain et al. (2016) #############################
############################################################

# PMF (Eq. 4)
dtdl <- function(x, p, beta){
  ((1-p)^2 * (1 + beta*x) * p^x) / (1 + p*(beta-1))
}

# Survival function (from paper)
Stdl <- function(x, p, beta){
  (p^x * ((1-p)*(1 + beta*x) + p*beta)) / (1 + p*(beta-1))
}

# CDF
Ftdl <- function(x, p, beta){
  1 - Stdl(x+1, p, beta)
}

# Log-likelihood
loglik.tdl <- function(par){
  p    <- par[1]
  beta <- par[2]
  if(p<=0 || p>=1 || beta<0) return(-Inf)
  sum(log(dtdl(x, p, beta)))
}

# MLE
ml.tdl <- maxLik(loglik.tdl, start=c(0.92, 227), method="NR")

p.tdl    <- ml.tdl$estimate[1]
beta.tdl <- ml.tdl$estimate[2]
se.tdl   <- c(sqrt(vcov(ml.tdl)[1,1]),sqrt(vcov(ml.tdl)[2,2]))

# CDF for KS
Fx.tdl <- function(x){
  Ftdl(x, p.tdl, beta.tdl)
}

# Goodness-of-fit
ks.tdl  <- ks.test(x, "Fx.tdl")
AIC.tdl <- -2*logLik(ml.tdl) + 2*2
BIC.tdl <- -2*logLik(ml.tdl) + 2*log(n)

############################################################
######## Poisson–XGamma (Altun et al. 2022) ################
############################################################

# PMF (Eq. 2)
dpxg <- function(x, theta){
  (theta^2 * (2*(1+theta)^2 + theta*(x+2)*(x+1))) /
    (2*(1+theta)^(x+4))
}

# CDF (from paper)
Fpxg <- function(x, theta){
  1 - (x^2*theta^2 + 5*x*theta^2 + 2*x*theta +
       2*theta^3 + 10*theta^2 + 8*theta + 2) /
      (2*(1+theta)^(x+4))
}

# Log-likelihood
loglik.pxg <- function(par){
  theta <- par[1]
  if(theta <= 0) return(-Inf)
  sum(log(dpxg(x, theta)))
}

# MLE
ml.pxg <- maxLik(loglik.pxg, start=c(1), method="NR")

theta.pxg <- ml.pxg$estimate
se.pxg    <- sqrt(diag(vcov(ml.pxg)))

# CDF for KS
Fx.pxg <- function(x){
  Fpxg(x, theta.pxg)
}

# Goodness-of-fit
ks.pxg  <- ks.test(x, "Fx.pxg")
AIC.pxg <- -2*logLik(ml.pxg) + 2
BIC.pxg <- -2*logLik(ml.pxg) + log(n)

############################################################
############ DsGLi (El-Morshedy et al. 2021) ################
############################################################
# PMF (Eq. 6)
ddsgli <- function(x, eta, alpha){
  if(any(eta <= 0 | eta >= 1 | alpha <= 0)) return(0)
  
  (eta^x / (1 - log(eta))) *
    (1 - eta
     - log(eta) * (1 + x - eta*(x + 2))
     + (1 - (1/alpha)) * (log(eta))^2 * (x - (x + 1)*eta))
}

# CDF (Eq. 7 – closed form)
Fdsgli <- function(x, eta, alpha){
  1 - (
    (
      ((1 - (x+1)*log(eta)) *
        (alpha - alpha*log(eta) + log(eta)) - log(eta)
      ) * eta^(x + 1)
    ) / (alpha * (1 - log(eta)))
  )
}

# Log-likelihood
loglik.dsgli <- function(par){
  eta   <- par[1]
  alpha <- par[2]
  if(eta <= 0 || eta >= 1 || alpha <= 0) return(-Inf)
  sum(log(ddsgli(x, eta, alpha)))
}

# MLE
ml.dsgli <- maxLik(loglik.dsgli, start = c(0.9, 5), method = "NR")

eta.dsgli   <- ml.dsgli$estimate[1]
alpha.dsgli <- ml.dsgli$estimate[2]
se.dsgli    <- sqrt(diag(vcov(ml.dsgli)))

# CDF for KS
Fx.dsgli <- function(x){
  Fdsgli(x, eta.dsgli, alpha.dsgli)
}

# GOF
ks.dsgli  <- ks.test(x, "Fx.dsgli")
AIC.dsgli <- -2 * logLik(ml.dsgli) + 2 * 2
BIC.dsgli <- -2 * logLik(ml.dsgli) + 2 * log(n)


############################################################
######## Discrete Power–Ailamujia (DsPA, 2022) ##############
############################################################

# PMF  (Eq. 6, page 3)
ddspa <- function(x, lambda, beta){
  if(any(lambda <= 0 | lambda >= 1 | beta <= 0)) return(0)
  
  lambda^(x^beta) * (1 - x^beta * log(lambda)) -
    lambda^((x + 1)^beta) * (1 - (x + 1)^beta * log(lambda))
}

# CDF  (Eq. 5, page 3 – closed form)
Fdspa <- function(x, lambda, beta){
  1 - lambda^((x + 1)^beta) * (1 - (x + 1)^beta * log(lambda))
}

# Log-likelihood
loglik.dspa <- function(par){
  lambda <- par[1]
  beta   <- par[2]
  if(lambda <= 0 || lambda >= 1 || beta <= 0) return(-Inf)
  
  sum(log(ddspa(x, lambda, beta)))
}

# MLE (as used in the paper)
ml.dspa <- maxLik(
  loglik.dspa,
  start  = c(0.5, 1),
  method = "NR"
)

lambda.dspa <- ml.dspa$estimate[1]
beta.dspa   <- ml.dspa$estimate[2]
se.dspa     <- sqrt(diag(vcov(ml.dspa)))

# CDF for KS test (closed form)
Fx.dspa <- function(x){
  Fdspa(x, lambda.dspa, beta.dspa)
}

# Goodness-of-fit
ks.dspa  <- ks.test(x, "Fx.dspa")
AIC.dspa <- -2 * logLik(ml.dspa) + 2 * 2
BIC.dspa <- -2 * logLik(ml.dspa) + 2 * log(n)

############################################################
######## Uniform–Geometric Distribution (UG) ###############
######## Akdoğan et al. (2016) #############################
############################################################

library(VGAM)   # for LerchPhi
library(maxLik)

# PMF (Eq. 2, page 4)
dug <- function(x, p){
  if(any(p <= 0 | p >= 1)) return(0)
  p * (1 - p)^(x - 1) * lerch(1 - p, 1, x)
}

# CDF (Eq. 6, page 6 – closed form)
Fug <- function(x, p){
  1 - p * (1 - p)^x * (1/p - x * lerch(1 - p, 1, x + 1))
}

# Log-likelihood
loglik.ug <- function(par){
  p <- par[1]
  if(p <= 0 || p >= 1) return(-Inf)
  sum(log(dug(x, p)))
}

# MLE (Newton–Raphson, as in paper)
ml.ug <- maxLik(
  loglik.ug,
  start  = 0.5,
  method = "NR"
)

p.ug  <- ml.ug$estimate
se.ug <- sqrt(vcov(ml.ug))

# CDF for KS test
Fx.ug <- function(x){
  Fug(x, p.ug)
}

# GOF
ks.ug  <- ks.test(x, "Fx.ug")
AIC.ug <- -2 * logLik(ml.ug) + 2 * 1
BIC.ug <- -2 * logLik(ml.ug) + log(n)

############################################################
###### Discrete Weibull–Geometric Distribution (DWG) #######
############################################################

library(maxLik)

# PMF  (Eq. 10 in the paper)
ddwg <- function(x, p, rho, alpha){
  if(any(p <= 0 | p >= 1 | rho <= 0 | rho >= 1 | alpha <= 0)) return(0)
  
  (1 - p) * (rho^(x^alpha) - rho^((x + 1)^alpha)) /
    ((1 - p * rho^(x^alpha)) * (1 - p * rho^((x + 1)^alpha)))
}

# CDF (Eq. 12 in the paper – closed form)
Fdwg <- function(x, p, rho, alpha){
  1 - rho^((x + 1)^alpha) / (1 - p * rho^((x + 1)^alpha))
}

# Log-likelihood
loglik.dwg <- function(par){
  p     <- par[1]
  rho   <- par[2]
  alpha <- par[3]
  
  if(p <= 0 || p >= 1 || rho <= 0 || rho >= 1 || alpha <= 0) return(-Inf)
  
  sum(log(ddwg(x, p, rho, alpha)))
}

# MLE (Newton–Raphson, as in paper)
ml.dwg <- maxLik(
  loglik.dwg,
  start  = c(0.5, 0.9, 1),
  method = "NR"
)

# Estimates
p.dwg     <- ml.dwg$estimate[1]
rho.dwg   <- ml.dwg$estimate[2]
alpha.dwg <- ml.dwg$estimate[3]

# Standard errors
se.dwg <- sqrt(diag(vcov(ml.dwg)))

# CDF for KS test
Fx.dwg <- function(q){
  Fdwg(q, p.dwg, rho.dwg, alpha.dwg)
}

# Goodness-of-fit
ks.dwg  <- ks.test(x, "Fx.dwg")
AIC.dwg <- -2 * logLik(ml.dwg) + 2 * 3
BIC.dwg <- -2 * logLik(ml.dwg) + 3 * log(n)

############################################################
############ BNPWE Distribution (Al-Bossly & Eliwa, 2022) ###
############################################################

library(maxLik)

# PMF  (Eq. 5)
dbnpwe <- function(x, alpha, beta, theta){
  if(any(alpha <= 0 | beta <= 0 | beta >= 1 | theta <= 0)) return(0)

  alpha * (1 + theta) * beta^x / (alpha + beta + alpha*theta)^(x + 1)
}

# CDF (Eq. 6 – closed form)
Fbnpwe <- function(x, alpha, beta, theta){
  1 - (beta / (alpha + beta + alpha*theta))^(x + 1)
}

# Log-likelihood
loglik.bnpwe <- function(par){
  alpha <- par[1]
  beta  <- par[2]
  theta <- par[3]

  if(alpha <= 0 || beta <= 0 || beta >= 1 || theta <= 0) return(-Inf)

  sum(log(dbnpwe(x, alpha, beta, theta)))
}

# MLE (Newton–Raphson)
ml.bnpwe <- maxLik(
  loglik.bnpwe,
  start  = c(0.5, 0.5, 0.5),
  method = "NR"
)

# Estimates
alpha.bnpwe <- ml.bnpwe$estimate[1]
beta.bnpwe  <- ml.bnpwe$estimate[2]
theta.bnpwe <- ml.bnpwe$estimate[3]

# Standard errors
se.bnpwe <- sqrt(diag(vcov(ml.bnpwe)))

# CDF for KS test
Fx.bnpwe <- function(q){
  Fbnpwe(q, alpha.bnpwe, beta.bnpwe, theta.bnpwe)
}

# GOF
ks.bnpwe  <- ks.test(x, "Fx.bnpwe")
AIC.bnpwe <- -2 * logLik(ml.bnpwe) + 2 * 3
BIC.bnpwe <- -2 * logLik(ml.bnpwe) + 3 * log(n)

############################################################
######## Discrete Teissier Distribution (DT) ################
############################################################

library(maxLik)

# PMF  (Eq. 5 in the paper)
ddt <- function(x, theta){
  if(any(theta <= 1)) return(0)
  exp(1) * theta^x * (exp(-theta^x) - theta * exp(-theta^(x + 1)))
}

# CDF  (Eq. 6 in the paper – closed form)
Fdt <- function(x, theta){
  1 - theta^(x + 1) * exp(1 - theta^(x + 1))
}

# Log-likelihood (log-scale, NaN-safe)
loglik.dt <- function(par){
  theta <- par[1]
  if(theta <= 1) return(-Inf)

  pmf <- ddt(x, theta)
  pmf[pmf <= 0] <- 1e-12  # numerical safety

  sum(log(pmf))
}

# MLE (stable optimizer)
ml.dt <- maxLik(
  loglik.dt,
  start  = c(1.2),
  method = "BFGS"
)

# Estimates
theta.dt <- ml.dt$estimate
se.dt    <- sqrt(diag(vcov(ml.dt)))

# CDF for KS test
Fx.dt <- function(q){
  Fdt(q, theta.dt)
}

# Goodness-of-fit
ks.dt  <- ks.test(x, "Fx.dt")
AIC.dt <- -2 * logLik(ml.dt) + 2 * 1
BIC.dt <- -2 * logLik(ml.dt) + log(n)

############################################################
##### Binomial Poisson–Ailamujia Distribution (BPA) #########
############################################################

library(maxLik)

# PMF  (Eq. 4 in the paper)
dBPA <- function(x, alpha, p){
  if(any(alpha <= 0) || any(p <= 0)) return(0)
  (alpha^2 * p^x * (1 + x)) / ((alpha + p)^(x + 2))
}

# CDF  (Eq. 2 in the paper – closed form)
FBPA <- function(x, alpha, p){
  1 - (p^(x + 1) * (2*alpha + alpha*x + p)) / ((alpha + p)^(x + 2))
}

# Log-likelihood (log-scale, NaN-safe)
loglik.BPA <- function(par){
  alpha <- par[1]
  p     <- par[2]
  
  if(alpha <= 0 || p <= 0) return(-Inf)
  
  pmf <- dBPA(x, alpha, p)
  pmf[pmf <= 0] <- 1e-12  # numerical safety
  
  sum(log(pmf))
}

# MLE (stable optimizer)
ml.BPA <- maxLik(
  loglik.BPA,
  start  = c(alpha = 1, p = 0.5),
  method = "BFGS"
)

# Estimates
alpha.BPA <- ml.BPA$estimate[1]
p.BPA     <- ml.BPA$estimate[2]
se.BPA    <- sqrt(diag(vcov(ml.BPA)))

# CDF for KS test
Fx.BPA <- function(q){
  FBPA(q, alpha.BPA, p.BPA)
}

# Goodness-of-fit
ks.BPA  <- ks.test(x, "Fx.BPA")
AIC.BPA <- -2 * logLik(ml.BPA) + 2 * 2
BIC.BPA <- -2 * logLik(ml.BPA) + log(n) * 2

############################################################
######## Discrete Moment Exponential (DME) Model ############
############################################################

library(maxLik)

# PMF (Eq. 3 in paper)
dDME <- function(x, theta){
  if(any(theta <= 0 | theta >= 1)) return(0)
  (1 - theta)^2 * x * theta^x
}

# CDF (Eq. 4 in paper)
FDME <- function(x, theta){
  1 - theta^x * (1 + x - theta * x)
}

# Log-likelihood (NaN-safe)
loglik.DME <- function(par){
  theta <- par[1]
  if(theta <= 0 || theta >= 1) return(-Inf)

  pmf <- dDME(x, theta)
  pmf[pmf <= 0] <- 1e-12

  sum(log(pmf))
}

# MLE (closed form is available, but we keep maxLik for uniformity)
ml.DME <- maxLik(
  loglik.DME,
  start  = c(0.5),
  method = "BFGS"
)

# Estimates
theta.DME <- ml.DME$estimate
se.DME    <- sqrt(diag(vcov(ml.DME)))

# CDF for KS test
Fx.DME <- function(q){
  FDME(q, theta.DME)
}

# Goodness-of-fit
ks.DME  <- ks.test(x, "Fx.DME")
AIC.DME <- -2 * logLik(ml.DME) + 2 * 1
BIC.DME <- -2 * logLik(ml.DME) + log(n) * 1

############################################################
###### Uniform Poisson–Ailamujia Distribution (UPA) #########
############################################################

library(maxLik)

# PMF  (Eq. 4)
dUPA <- function(x, alpha){
  if(any(alpha <= 0)) return(0)
  (2 * alpha) / (1 + 2 * alpha)^(x + 1)
}

# CDF  (Eq. 6)
FUPA <- function(x, alpha){
  1 - 1 / (1 + 2 * alpha)^(x + 1)
}

# Log-likelihood (safe)
loglik.UPA <- function(par){
  alpha <- par[1]
  if(alpha <= 0) return(-Inf)
  
  pmf <- dUPA(x, alpha)
  pmf[pmf <= 0] <- 1e-12
  sum(log(pmf))
}

# MLE (BFGS works perfectly here)
ml.UPA <- maxLik(
  loglik.UPA,
  start  = c(alpha = 1),
  method = "BFGS"
)

# Estimates
alpha.UPA <- ml.UPA$estimate
se.UPA    <- sqrt(diag(vcov(ml.UPA)))

# CDF for KS test
Fx.UPA <- function(q){
  FUPA(q, alpha.UPA)
}

# Goodness-of-fit
ks.UPA  <- ks.test(x, "Fx.UPA")
AIC.UPA <- -2 * logLik(ml.UPA) + 2
BIC.UPA <- -2 * logLik(ml.UPA) + log(n)

############################################################
###### Discrete Inverted Nadarajah–Haghighi (DINH) ##########
############################################################

library(maxLik)

# PMF  (Eq. 3 in paper)
dinh <- function(x, alpha, lambda){
  if(any(alpha <= 0 | lambda <= 0)) return(0)

  exp(1) * (
    exp(-(1 + lambda/(x+1))^alpha) -
      exp(-(1 + lambda/x)^alpha)
  )
}

# CDF (Eq. 4 in paper – closed form)
Finh <- function(x, alpha, lambda){
  exp(1 - (1 + lambda/(x+1))^alpha)
}

# Log-likelihood (numerically safe)
loglik.dinh <- function(par){
  alpha  <- par[1]
  lambda <- par[2]
  if(alpha <= 0 || lambda <= 0) return(-Inf)

  pmf <- dinh(x, alpha, lambda)
  pmf[pmf <= 0] <- 1e-12
  sum(log(pmf))
}

# MLE (Newton–Raphson as in paper)
ml.dinh <- maxLik(
  loglik.dinh,
  start  = c(alpha = 2, lambda = 4),
  method = "NM"
)

# Estimates
alpha.dinh  <- ml.dinh$estimate[1]
lambda.dinh <- ml.dinh$estimate[2]

# Standard errors
se.dinh <- sqrt(diag(vcov(ml.dinh)))

# CDF for KS test
Fx.dinh <- function(q){
  Finh(q, alpha.dinh, lambda.dinh)
}

# Goodness-of-fit
ks.dinh  <- ks.test(x, "Fx.dinh")
AIC.dinh <- -2 * logLik(ml.dinh) + 2 * 2
BIC.dinh <- -2 * logLik(ml.dinh) + 2 * log(n)

#########################################################
-Log Likelihood for all Fitted Models
########################################################
NegLogLik <- c(
  -logLik(ml.bdme),
  -logLik(ml.pois),
  -logLik(ml.bin),
  -logLik(ml.sgeom),
  -logLik(ml.nb),
  -logLik(ml.dp),
  -logLik(ml.dw),
  -logLik(ml.gp),
  -logLik(ml.dlog),
  -logLik(ml.dburr),
  -logLik(ml.dg),
  -logLik(ml.dlin2),
  -logLik(ml.edl),
  -logLik(ml.eg),
  -logLik(ml.tdl),
  -logLik(ml.pxg),
  -logLik(ml.dsgli),
  -logLik(ml.dspa),
  -logLik(ml.ug),
  -logLik(ml.dinh),
  -logLik(ml.dwg),
  -logLik(ml.bnpwe),
  -logLik(ml.dt),
  -logLik(ml.BPA),
  -logLik(ml.DME),
  -logLik(ml.UPA)
)
############################################################
############################################################
######################## Results Table #####################
############################################################
results <- data.frame(
  Model = c("BDME","Poisson","Binomial","Shifted Geometric",
            "Negative Binomial","Discrete Perks","Discrete Weibull",
            "Generalized Poisson","Discrete Logistic","Discrete Burr",
            "Discrete Gamma","Discrete Lindley II",
            "Exponentiated Discrete Lindley","Exponentiated Geometric",
            "Two Parameter Discrete Lindley","Poisson–XGamma",
            "DsGLi","DsPA","Uniform–Geometric","DINH",
            "Discrete Weibull–Geometric","BNPWE","Discrete Teissier",
            "Binomial Poisson–Ailamujia","DME",
            "Uniform Poisson–Ailamujia"),

  MLE_SE = c(
    paste0(round(p.bdme,4)," (",round(se.bdme,4),")"),
    paste0(round(lam.hat,4)," (",round(se.pois,4),")"),
    paste0(round(p.bin,4)," (",round(se.bin,4),")"),
    paste0(round(p.sg,4)," (",round(se.sg,4),")"),
    paste0("r=",round(r.nb,4),", p=",round(p.nb,4),
           " (",paste(round(se.nb,4),collapse=", "),")"),
    paste0("α=",round(alpha.dp,4),", β=",round(beta.dp,4),
           " (",paste(round(se.dp,4),collapse=", "),")"),
    paste0("β=",round(beta.dw,4),", p=",round(p.dw,4),
           " (",paste(round(se.dw,4),collapse=", "),")"),
    paste0("λ1=",round(lambda1.gp,4),", λ2=",round(lambda2.gp,4),
           " (",paste(round(se.gp,4),collapse=", "),")"),
    paste0("p=",round(p.dlog,4),", μ=",round(mu.dlog,4),
           " (",paste(round(se.dlog,4),collapse=", "),")"),
    paste0("θ=",round(theta.db,4),", α=",round(alpha.db,4),
           " (",paste(round(se.dburr,4),collapse=", "),")"),
    paste0("shape=",round(shape.dg,4),", scale=",round(scale.dg,4),
           " (",paste(round(se.dg,4),collapse=", "),")"),
    paste0("θ=",round(theta.dlin2,4)," (",round(se.dlin2,4),")"),
    paste0("a=",round(a.edl,4),", b=",round(b.edl,4),
           " (",paste(round(se.edl,4),collapse=", "),")"),
    paste0("q=",round(q.eg,4),", α=",round(alpha.eg,4),
           " (",paste(round(se.eg,4),collapse=", "),")"),
    paste0("p=",round(p.tdl,4),", β=",round(beta.tdl,4),
           " (",paste(round(se.tdl,4),collapse=", "),")"),
    paste0("θ=",round(theta.pxg,4)," (",round(se.pxg,4),")"),
    paste0("η=",round(eta.dsgli,4),", α=",round(alpha.dsgli,4),
           " (",paste(round(se.dsgli,4),collapse=", "),")"),
    paste0("λ=",round(lambda.dspa,4),", β=",round(beta.dspa,4),
           " (",paste(round(se.dspa,4),collapse=", "),")"),
    paste0("p=",round(p.ug,4)," (",round(se.ug,4),")"),
paste0("α=",round(alpha.dinh,4),", λ=",round(lambda.dinh,4),
           " (",paste(round(se.dinh,4),collapse=", "),")"),
    paste0("p=",round(p.dwg,4),", ρ=",round(rho.dwg,4),", α=",round(alpha.dwg,4),
           " (",paste(round(se.dwg,4),collapse=", "),")"),
    paste0("α=",round(alpha.bnpwe,4),", β=",round(beta.bnpwe,4),
           ", θ=",round(theta.bnpwe,4),
           " (",paste(round(se.bnpwe,4),collapse=", "),")"),
    paste0("θ=",round(theta.dt,4)," (",round(se.dt,4),")"),
    paste0("α=",round(alpha.BPA,4),", p=",round(p.BPA,4),
           " (",paste(round(se.BPA,4),collapse=", "),")"),
    paste0("θ=",round(theta.DME,4)," (",round(se.DME,4),")"),
    paste0("α=",round(alpha.UPA,4)," (",round(se.UPA,4),")")
  ),

NegLogLik = as.numeric(NegLogLik),

  AIC = c(AIC.bdme, AIC.pois, AIC.bin, AIC.sgeom,
          AIC.nb, AIC.dp, AIC.dw, AIC.gp, AIC.dlog,
          AIC.dburr, AIC.dg, AIC.dlin2,
          AIC.edl, AIC.eg, AIC.tdl, AIC.pxg,
          AIC.dsgli, AIC.dspa, AIC.ug,  AIC.dinh, AIC.dwg,
          AIC.bnpwe, AIC.dt, AIC.BPA, AIC.DME,
          AIC.UPA),

  BIC = c(BIC.bdme, BIC.pois, BIC.bin, BIC.sgeom,
          BIC.nb, BIC.dp, BIC.dw, BIC.gp, BIC.dlog,
          BIC.dburr, BIC.dg, BIC.dlin2,
          BIC.edl, BIC.eg, BIC.tdl, BIC.pxg,
          BIC.dsgli, BIC.dspa, BIC.ug, BIC.dinh,BIC.dwg,
          BIC.bnpwe, BIC.dt, BIC.BPA, BIC.DME,
          BIC.UPA),

  KS = c(ks.bdme$statistic, ks.pois$statistic, ks.bin$statistic,
         ks.sgeom$statistic, ks.nb$statistic, ks.dp$statistic,
         ks.dw$statistic, ks.gp$statistic, ks.dlog$statistic,
         ks.dburr$statistic, ks.dg$statistic,
         ks.dlin2$statistic, ks.edl$statistic, ks.eg$statistic,
         ks.tdl$statistic, ks.pxg$statistic, ks.dsgli$statistic,
         ks.dspa$statistic, ks.ug$statistic, ks.dinh$statistic, ks.dwg$statistic,
         ks.bnpwe$statistic, ks.dt$statistic, ks.BPA$statistic,
         ks.DME$statistic, ks.UPA$statistic),

  KS_p = c(ks.bdme$p.value, ks.pois$p.value, ks.bin$p.value,
           ks.sgeom$p.value, ks.nb$p.value, ks.dp$p.value,
           ks.dw$p.value, ks.gp$p.value, ks.dlog$p.value,
           ks.dburr$p.value, ks.dg$p.value,
           ks.dlin2$p.value, ks.edl$p.value, ks.eg$p.value,
           ks.tdl$p.value, ks.pxg$p.value, ks.dsgli$p.value,
           ks.dspa$p.value, ks.ug$p.value, ks.dinh$p.value, ks.dwg$p.value,
           ks.bnpwe$p.value, ks.dt$p.value, ks.BPA$p.value,
           ks.DME$p.value, ks.UPA$p.value)
)

print(results)
