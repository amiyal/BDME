rm(list = ls(all = TRUE))
library(maxLik)
library(plyr)
set.seed(2025)

p0 <- 0.75
NN <- 1000
n_iter <- 110000
burn_in <- 10000
init_p <- p0
proposal_sd <- 0.05
a_prior <- 9
b_prior <- 3
target_accept <- 0.3

n_vals <- c(10, 20, 30, 50, 100, 200)

AE.mle <- MSE.mle <- Bias.mle <- numeric(length(n_vals))
AE.mo  <- MSE.mo  <- Bias.mo  <- numeric(length(n_vals))
AE.mps <- MSE.mps <- Bias.mps <- numeric(length(n_vals))
AE.po <- MSE.po <- Bias.po <- numeric(length(n_vals))
AE.be.ip <- MSE.be.ip <- Bias.be.ip <- numeric(length(n_vals))
AE.be.nip <- MSE.be.nip <- Bias.be.nip <- numeric(length(n_vals))


# Function to compute metrics
compute_metrics <- function(estimates, true_val) {
  c(
    AE = mean(estimates, na.rm = TRUE),
    MSE = mean((estimates - true_val)^2, na.rm = TRUE),
    Bias = abs(mean(estimates - true_val, na.rm = TRUE))
  )
}


# Define PMF and CDF for BNDME
dBNDME <- function(x, p) {
  ((1 - p)^2) * (p^(2 * x - 1)) * (x + p * (1 - p)) / ((1 - p * (1 - p))^(x + 2))
}

pBNDME <- function(x, p) {
  1 - (p^(2 * x + 1)) * (1 + (1 - p) * (x - p^2)) / ((1 - p * (1 - p))^(x + 2))
}

log_likelihood <- function(data, p) {
  if (p <= 0 || p >= 1) return(-Inf)
  sum(log(sapply(data, function(x) dBNDME(x, p))))
}

log_prior <- function(p, a , b ) {
  if (p <= 0 || p >= 1) return(-Inf)
  dbeta(p, a, b, log = TRUE)
}

##Bayes

# ==== MH Sampler ====
log_prior <- function(p, a , b ) {
  if (p <= 0 || p >= 1) return(-Inf)
  dbeta(p, a, b, log = TRUE)
}

log_posterior <- function(p, data, a , b) {
  if (p <= 0 || p >= 1) return(-Inf)
  sum(log(sapply(data, function(x) dBNDME(x, p)))) + log_prior(p, a, b)
}

mh_sampler_adaptive <- function(data, n_iter, init, proposal_sd, a, b, target_accept = 0.3) {
  samples <- numeric(n_iter)
  samples[1] <- init
  accept <- 0
  sd_curr <- proposal_sd
  
  for (i in 2:n_iter) {
    current <- samples[i - 1]
    proposal <- rnorm(1, mean = current, sd = sd_curr)
    if (proposal <= 0 || proposal >= 1) {
      samples[i] <- current
      next
    }
    
    log_alpha <- log_posterior(proposal, data, a, b) - log_posterior(current, data, a, b)
    alpha <- exp(log_alpha)
    
    if (runif(1) < alpha) {
      samples[i] <- proposal
      accept <- accept + 1
    } else {
      samples[i] <- current
    }
    
    if (i %% 100 == 0 && i <= n_iter / 2) {
      acc_rate <- accept / i
      if (acc_rate > target_accept + 0.05) {
        sd_curr <- sd_curr * 1.1
      } else if (acc_rate < target_accept - 0.05) {
        sd_curr <- sd_curr * 0.9
      }
    }
  }
  
  return(samples)
}

for (k in seq_along(n_vals)) {
  n_k <- n_vals[k]
  p.hat.mle <- p.hat.mo <- p.hat.mps <- p.hat.po<-p.hat.be.ip<-p.hat.be.nip<- numeric(NN)
  
  for (ii in 1:NN) {
    # Simulate y from inverse CDF
    y <- numeric(n_k)
    for (i in 1:n_k) {
      u <- runif(1)
      q <- function(x) {
        (2 * x + 1) * log(p0) +
          log(1 + (1 - p0) * (x - p0^2)) -
          (x + 2) * log(1 - p0 * (1 - p0)) -
          log(1 - u)
      }
      y[i] <- tryCatch(uniroot(q, c(-0.5, 200), extendInt = "yes")$root, error = function(e) NA)
    }
    
    if (any(is.na(y))) next  # skip iteration if root-finding fails
    
    x <- ceiling(y)
    xbar <- mean(x)
    
    ## --- MLE Estimation ---
    loglik <- function(theta) {
      p <- theta[1]
      if (p <= 0 || p >= 1 || any(x + p * (1 - p) <= 0) || (1 - p * (1 - p)) <= 0) return(-Inf)
      2 * n_k * log(1 - p) +
        (2 * n_k * xbar - n_k) * log(p) +
        sum(log(x + p * (1 - p))) -
        (n_k * xbar + 2 * n_k) * log(1 - p * (1 - p))
    }
    
    ml <- tryCatch(maxLik(loglik, start = 0.05), error = function(e) NULL)
    if (!is.null(ml) && !any(is.na(ml$estimate))) {
      p.hat.mle[ii] <- ml$estimate
    } else {
      p.hat.mle[ii] <- NA
    }
    
    ## --- MOM Estimator ---
    p.hat.mo[ii] <- (-(xbar + 1) + sqrt(xbar^2 + 6 * xbar + 1)) / 2
    
    ## --- MPS Estimation ---
    x <- sort(x)
    fmps<- function(ppt){
      p<-ppt[1]
      D<-NULL
      D[1]<-pBNDME(x[1],p)
      D[n_k+1]<-1-pBNDME(x[n_k],p)
      for (i in 2:n_k){
        if(x[i]==x[i-1]){
          D[i]<- dBNDME(x[i],p)
          next
        } 
        D[i]<-pBNDME(x[i],p)-pBNDME(x[i-1],p)
      }
      aux<-1/(n_k+1)*sum(log(D))
      return(aux)
    }                 
    
    amps <- tryCatch(maxBFGS(fmps, start = c(p0))$estimate, error = function(e) NA)
    if (is.double(amps[1]) && amps[1] > 0 && amps[1] < 1) {
      p.hat.mps[ii] <- amps
    } else {
      p.hat.mps[ii] <- NA
    }
    
    
    # Proportion-of-zeros estimator
    estimate_p_pz <- function(x, n) {
      tab <- count(x)
      N0 <- ifelse(0 %in% tab$x, tab$freq[tab$x == 0], 0)
      f0 <- N0 / n
      if (f0 == 0 || is.na(f0)) return(NA)
      
      fp <- function(p) {
        (1 - p)^3 - ((1 - p * (1 - p))^2) * f0
      }
      
      res <- tryCatch({
        uniroot(fp, c(0.1, 0.99), extendInt = "yes")$root
      }, error = function(e) {
        NA
      })
      return(res)
    }
    
    p_pzo <- estimate_p_pz(x, n_k)
    p.hat.po[ii]<- p_pzo
    
    
    samp.ip <- mh_sampler_adaptive(x, n_iter, init_p, proposal_sd, a_prior, b_prior)
    posterior.ip <- samp.ip[seq(burn_in + 1, n_iter, by = 10)]
    p.hat.be.ip[ii] <- mean(posterior.ip)
    
    
    samp.nip <- mh_sampler_adaptive(x, n_iter, init_p, proposal_sd, a=1, b=1)
    posterior.nip <- samp.nip[seq(burn_in + 1, n_iter, by = 10)]
    p.hat.be.nip[ii] <- mean(posterior.nip)
    
    cat("  Iteration:", ii, "of", NN, "for sample size", n_k, "\n")
    
  } # End of NN replications
  
  # summary
  AE.mle[k] <- compute_metrics(p.hat.mle, p0)[1]
  MSE.mle[k] <- compute_metrics(p.hat.mle, p0)[2]
  Bias.mle[k] <- compute_metrics(p.hat.mle, p0)[3]
  
  AE.mo[k] <- compute_metrics(p.hat.mo, p0)[1]
  MSE.mo[k] <- compute_metrics(p.hat.mo, p0)[2]
  Bias.mo[k] <- compute_metrics(p.hat.mo, p0)[3]
  
  AE.mps[k] <- compute_metrics(p.hat.mps, p0)[1]
  MSE.mps[k] <- compute_metrics(p.hat.mps, p0)[2]
  Bias.mps[k] <- compute_metrics(p.hat.mps, p0)[3]
  
  AE.po[k] <- compute_metrics(p.hat.po, p0)[1]
  MSE.po[k] <- compute_metrics(p.hat.po, p0)[2]
  Bias.po[k] <- compute_metrics(p.hat.po, p0)[3]
  
  AE.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[1]
  MSE.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[2]
  Bias.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[3]
  
  AE.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[1]
  MSE.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[2]
  Bias.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[3]
  
  cat("Completed sample size:", n_k, "\n")
}

# Final output
result_df <- data.frame(
  Sample_Size = n_vals,
  AE_MLE = AE.mle, MSE_MLE = MSE.mle, Bias_MLE = Bias.mle,
  AE_MOM = AE.mo, MSE_MOM = MSE.mo, Bias_MOM = Bias.mo,
  AE_MPS = AE.mps, MSE_MPS = MSE.mps, Bias_MPS = Bias.mps,
  AE_PO = AE.po, MSE_PO = MSE.po, Bias_PO = Bias.po,
  AE_BE_IP = AE.be.ip, MSE_BE_IP = MSE.be.ip, Bias_BE_IP = Bias.be.ip,
  AE_BE_NIP = AE.be.nip, MSE_BE_NIP = MSE.be.nip, Bias_BE_NIP = Bias.be.nip
  
)

print(result_df)
write.csv(result_df, "C:/Users/Admin/Desktop/Prof.Bhupendra/BNDME_Estimation_Results1.csv", row.names = FALSE) #p=0.25
write.csv(result_df, "C:/Users/Admin/Desktop/Prof.Bhupendra/BNDME_Estimation_Results2.csv", row.names = FALSE)#p=0.5
write.csv(result_df, "C:/Users/Admin/Desktop/Prof.Bhupendra/BNDME_Estimation_Results3.csv", row.names = FALSE) #p=0.75

