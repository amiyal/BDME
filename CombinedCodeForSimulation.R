rm(list = ls(all = TRUE))
library(maxLik)
library(plyr)
set.seed(2024)


#p0 <- 0.10; a_prior <- 1; b_prior <- 9
#p0 <- 0.25; a_prior <- 3; b_prior <- 9
#p0 <- 0.5; a_prior <- 3; b_prior <- 3
#p0 <- 0.75; a_prior <- 9; b_prior <- 3
p0 <- 0.90; a_prior <- 9; b_prior <- 1

NN <- 1000
n_iter <- 110000
burn_in <- 10000
init_p <- p0
proposal_sd <- 0.05

target_accept <- 0.3
n_vals <- c(10, 20, 30, 50, 70, 100, 200, 500)
            
AE.mle      <- MSE.mle    <- AB.mle    <- Bias.mle    <- MRE.mle    <- numeric(length(n_vals))
AE.mo       <- MSE.mo     <- AB.mo     <- Bias.mo     <- MRE.mo     <- numeric(length(n_vals))
AE.mps      <- MSE.mps    <- AB.mps    <- Bias.mps    <- MRE.mps    <- numeric(length(n_vals))
AE.po       <- MSE.po     <- AB.po     <- Bias.po     <- MRE.po     <- numeric(length(n_vals))
AE.be.ip    <- MSE.be.ip  <- AB.be.ip  <- Bias.be.ip  <- MRE.be.ip  <- numeric(length(n_vals))
AE.be.nip   <- MSE.be.nip <- AB.be.nip <- Bias.be.nip <- MRE.be.nip <- numeric(length(n_vals))


# Function to compute metrics
compute_metrics <- function(estimates, true_val) {
  c(
    AE = mean(estimates, na.rm = TRUE),
    MSE = mean((estimates - true_val)^2, na.rm = TRUE),
    AB = mean(abs(estimates - true_val), na.rm = TRUE),
    Bias= mean(estimates - true_val, na.rm = TRUE),
    MRE= mean(estimates / true_val, na.rm = TRUE)
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
  sum(log(dBNDME(data, p)))
}

log_posterior <- function(p, data, a, b) {
  if (p <= 0 || p >= 1) return(-Inf)
  sum(log(dBNDME(data, p))) + dbeta(p, a, b, log = TRUE)
}

mh_sampler_adaptive <- function(data, n_iter, init, proposal_sd, a, b, target_accept = 0.3) {
  
  samples <- numeric(n_iter)
  samples[1] <- init
  accept <- 0
  sd_curr <- proposal_sd
  
  logpost_curr <- log_posterior(init, data, a, b)
  
  for (i in 2:n_iter) {
    
    proposal <- rnorm(1, samples[i - 1], sd_curr)
    
    if (proposal <= 0 || proposal >= 1) {
      samples[i] <- samples[i - 1]
      next
    }
    
    logpost_prop <- log_posterior(proposal, data, a, b)
    log_alpha <- logpost_prop - logpost_curr
    
    if (runif(1) < exp(log_alpha)) {
      samples[i] <- proposal
      logpost_curr <- logpost_prop
      accept <- accept + 1
    } else {
      samples[i] <- samples[i - 1]
    }
    
    if (i %% 100 == 0 && i <= n_iter / 2) {
      acc_rate <- accept / i
      sd_curr <- if (acc_rate > target_accept + 0.05) sd_curr * 1.1
      else if (acc_rate < target_accept - 0.05) sd_curr * 0.9
      else sd_curr
    }
  }
  samples
}

for (k in seq_along(n_vals)) {
  n_k <- n_vals[k]
  p.hat.mle <- p.hat.mo <- p.hat.mps <- p.hat.po <- p.hat.be.ip <- p.hat.be.nip <- numeric(NN)
  
  for (ii in 1:NN) {
    # Simulate y from inverse CDF
    y <- numeric(n_k)
    x <- numeric(n_k)
    for (i in 1:n_k) {
      u=runif(1,0,1)
      y[i] <- ceiling(uniroot(function(z) 1 - p0^z * (1 + z - p0*z) - u,c(0, 200), extendInt = "yes")$root)
      x[i]=rbinom(1,y[i],p0)
    }
    
    if (any(is.na(x))) next  # skip iteration if root-finding fails
    
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
    xs <- sort(x)
    D <- numeric(n_k + 1)
    fmps <- function(ppt) {
      
      p <- ppt[1]
      if (p <= 0 || p >= 1) return(-Inf)
      
      Fx <- pBNDME(xs, p)      # vectorized CDF
      fx <- dBNDME(xs, p)      # vectorized PDF
      
      D[1] <- Fx[1]
      D[n_k + 1] <- 1 - Fx[n_k]
      
      for (i in 2:n_k) {
        if (xs[i] == xs[i - 1]) {
          D[i] <- fx[i]
        } else {
          D[i] <- Fx[i] - Fx[i - 1]
        }
      }
      
      mean(log(D))
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
    
      cat("Iteration:", ii, "Sample size:", n_k, "\n")
    
  } # End of NN replications
  
  # summary
  AE.mle[k] <- compute_metrics(p.hat.mle, p0)[1]
  MSE.mle[k] <- compute_metrics(p.hat.mle, p0)[2]
  AB.mle[k] <- compute_metrics(p.hat.mle, p0)[3]
  Bias.mle[k] <- compute_metrics(p.hat.mle, p0)[4]
  MRE.mle[k] <- compute_metrics(p.hat.mle, p0)[5]
  
  
  AE.mo[k] <- compute_metrics(p.hat.mo, p0)[1]
  MSE.mo[k] <- compute_metrics(p.hat.mo, p0)[2]
  AB.mo[k] <- compute_metrics(p.hat.mo, p0)[3]
  Bias.mo[k] <- compute_metrics(p.hat.mo, p0)[4]
  MRE.mo[k] <- compute_metrics(p.hat.mo, p0)[5]
  
  AE.mps[k] <- compute_metrics(p.hat.mps, p0)[1]
  MSE.mps[k] <- compute_metrics(p.hat.mps, p0)[2]
  AB.mps[k] <- compute_metrics(p.hat.mps, p0)[3]
  Bias.mps[k] <- compute_metrics(p.hat.mps, p0)[4]
  MRE.mps[k] <- compute_metrics(p.hat.mps, p0)[5]
  
  
  AE.po[k] <- compute_metrics(p.hat.po, p0)[1]
  MSE.po[k] <- compute_metrics(p.hat.po, p0)[2]
  AB.po[k] <- compute_metrics(p.hat.po, p0)[3]
  Bias.po[k] <- compute_metrics(p.hat.po, p0)[4]
  MRE.po[k] <- compute_metrics(p.hat.po, p0)[5]
  
  AE.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[1]
  MSE.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[2]
  AB.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[3]
  Bias.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[4]
  MRE.be.ip[k] <- compute_metrics(p.hat.be.ip, p0)[5]
  
  AE.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[1]
  MSE.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[2]
  AB.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[3]
  Bias.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[4]
  MRE.be.nip[k] <- compute_metrics(p.hat.be.nip, p0)[5]
  
  cat("Completed sample size:", n_k, "\n")
}

# Final output
result_df <- data.frame(
  Sample_Size = n_vals,
  AE_MLE = AE.mle, MSE_MLE = MSE.mle, AB_MLE= AB.mle, Bias_MLE = Bias.mle, MRE_MLE= MRE.mle,
  AE_MOM = AE.mo, MSE_MOM = MSE.mo, AB_MOM= AB.mo, Bias_MOM = Bias.mo, MRE_MOM= MRE.mo,
  AE_MPS = AE.mps, MSE_MPS = MSE.mps, AB_MPS= AB.mps, Bias_MPS = Bias.mps, MRE_MPS= MRE.mps,
  AE_PO = AE.po, MSE_PO = MSE.po, AB_PO= AB.po,  Bias_PO = Bias.po, MRE_PO= MRE.po,
  AE_BE_IP = AE.be.ip,  MSE_BE_IP = MSE.be.ip, AB_BE_IP= AB.be.ip,  Bias_BE_IP = Bias.be.ip, MRE_BE_IP= MRE.be.ip,
  AE_BE_NIP = AE.be.nip, MSE_BE_NIP = MSE.be.nip, AB_BE_NIP= AB.be.nip,  Bias_BE_NIP = Bias.be.nip, MRE_BE_NIP= MRE.be.nip
  
)

print(result_df)
